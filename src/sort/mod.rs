use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use anyhow::{Context, Result};
use bstr::BString;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as AlignmentWrite;
use noodles_sam::header::record::value::Map;
use noodles_sam::header::record::value::map::header::{Version, tag::SORT_ORDER};
use rayon::prelude::*;
use tempfile::TempDir;

use crate::cli::SortArgs;
use crate::io::{
    BamReader, BamWriter, OpenOptions, PgInfo, WriteOptions, append_pg, install_thread_pool,
    is_record_level_error, resolve_memory_hint,
};
use crate::markdup::{MarkdupOptions, mark_duplicates_in_memory};
use crate::types::OutputFormat;

/// Presizes the chunk vector only; spilling is decided by measured bytes.
const ESTIMATED_BYTES_PER_RECORD: usize = 512;

/// Approximate live heap bytes of one decoded record. Budgeting in records
/// instead would overshoot by two orders of magnitude on long reads, where a
/// record costs ~200 kB rather than a few hundred bytes.
fn approx_heap_size(rec: &RecordBuf) -> usize {
    /// Struct, `Vec` headers and allocator round-up.
    const FIXED: usize = 160;
    /// Undercounts long tag values.
    const PER_TAG: usize = 24;

    FIXED
        + rec.name().map(|n| n.len()).unwrap_or(0)
        + rec.sequence().len()
        + rec.quality_scores().as_ref().len()
        + rec.cigar().as_ref().len() * 4
        + rec.data().len() * PER_TAG
}

/// Minimum sort buffer when `-m auto` resolves below this. Picked so the
/// 0.5 M-record short batches still get a useful in-memory chunk.
const AUTO_MIN_MEMORY_BYTES: usize = 512 * 1024 * 1024; // 512 MB

/// `auto` allocates 3/4 of total RAM. Matches the user's request — on a
/// 32 GB box we end up with 24 GB sort budget, comfortably hosting
/// chr20 30× WGS in a single chunk.
const AUTO_MEMORY_NUM: u32 = 3;
const AUTO_MEMORY_DEN: u32 = 4;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SortKind {
    Coordinate,
    QueryName,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum SortKey {
    /// Coordinate order, matching `samtools sort`. htslib's key packs
    /// `(tid, pos, reverse-strand-bit)` — i.e. at one position forward reads sort
    /// before reverse reads — and then a *stable* sort keeps equal-key records in
    /// input order. (We used to fold the whole flag + qname into the key, which
    /// both diverged from samtools and allocated a qname per record for nothing.)
    Coord { tid: u32, pos: u32, rev: u8 },
    /// Name order, matching `samtools sort -n`: natural/`strnum` qname comparison,
    /// then the READ1/READ2 bits (`flag & 0xc0`); a stable sort handles the rest.
    Name { qname: Vec<u8>, read12: u8 },
}

impl Ord for SortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (
                SortKey::Coord {
                    tid: a_tid,
                    pos: a_pos,
                    rev: a_rev,
                },
                SortKey::Coord {
                    tid: b_tid,
                    pos: b_pos,
                    rev: b_rev,
                },
            ) => (a_tid, a_pos, a_rev).cmp(&(b_tid, b_pos, b_rev)),
            (
                SortKey::Name {
                    qname: a_q,
                    read12: a_r,
                },
                SortKey::Name {
                    qname: b_q,
                    read12: b_r,
                },
            ) => strnum_cmp(a_q, b_q).then_with(|| a_r.cmp(b_r)),
            _ => unreachable!("mixed sort kinds"),
        }
    }
}

impl PartialOrd for SortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Natural "string + number" comparison, a faithful port of htslib's
/// `strnum_cmp` (used by `samtools sort -n`): maximal digit runs compare as
/// numbers (leading zeros ignored, then more-leading-zeros sorts greater),
/// everything else byte-by-byte. This makes `r2` sort before `r10` exactly as
/// samtools does.
pub(crate) fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (na, nb) = (a.len(), b.len());
    let (mut ia, mut ib) = (0usize, 0usize);
    while ia < na && ib < nb {
        if a[ia].is_ascii_digit() && b[ib].is_ascii_digit() {
            while ia < na && a[ia] == b'0' {
                ia += 1;
            }
            while ib < nb && b[ib] == b'0' {
                ib += 1;
            }
            // Walk the common run of identical digits.
            while ia < na
                && ib < nb
                && a[ia].is_ascii_digit()
                && b[ib].is_ascii_digit()
                && a[ia] == b[ib]
            {
                ia += 1;
                ib += 1;
            }
            let da = ia < na && a[ia].is_ascii_digit();
            let db = ib < nb && b[ib].is_ascii_digit();
            if da && db {
                // First differing digit: whichever has the longer remaining digit
                // run is the larger number.
                let mut k = 0usize;
                loop {
                    let pa = ia + k < na && a[ia + k].is_ascii_digit();
                    let pb = ib + k < nb && b[ib + k].is_ascii_digit();
                    if pa && pb {
                        k += 1;
                        continue;
                    }
                    return if pa {
                        Ordering::Greater
                    } else if pb {
                        Ordering::Less
                    } else {
                        a[ia].cmp(&b[ib])
                    };
                }
            } else if da {
                return Ordering::Greater;
            } else if db {
                return Ordering::Less;
            } else if ia != ib {
                // Equal value but a different number of leading zeros consumed.
                return ia.cmp(&ib);
            }
            // Equal numbers, same offset — fall through and keep comparing.
        } else if a[ia] != b[ib] {
            return a[ia].cmp(&b[ib]);
        } else {
            ia += 1;
            ib += 1;
        }
    }
    if ia < na {
        Ordering::Greater
    } else if ib < nb {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

/// Sort `records` by the requested key, in parallel.
///
/// Sorting compact `(key, index)` pairs and permuting once keeps the ~200-byte
/// `RecordBuf`s out of the `n log n` swaps, and avoids the qname allocation a
/// key-rebuilding comparator would make per comparison.
///
/// The record index is the last component of every key, so ties keep input order
/// as samtools does and an unstable sort suffices.
fn sort_records_by_kind(records: &mut Vec<RecordBuf>, kind: SortKind) {
    match kind {
        SortKind::Coordinate => {
            let mut keys: Vec<(u32, u32, u8, u32)> = records
                .par_iter()
                .enumerate()
                .map(|(i, r)| {
                    let (tid, pos, rev) = coord_key(r);
                    (tid, pos, rev, i as u32)
                })
                .collect();
            keys.par_sort_unstable();
            let order: Vec<u32> = keys.iter().map(|k| k.3).collect();
            apply_permutation(records, order);
        }
        SortKind::QueryName => {
            // One arena, so a key is a slice range rather than an allocation.
            let mut arena: Vec<u8> = Vec::with_capacity(records.len() * 32);
            let mut keys: Vec<NameKey> = Vec::with_capacity(records.len());
            for (i, r) in records.iter().enumerate() {
                let name: &[u8] = r.name().map(|n| n.as_ref()).unwrap_or_default();
                let start = arena.len() as u32;
                arena.extend_from_slice(name);
                keys.push(NameKey {
                    start,
                    len: name.len() as u32,
                    read12: (u16::from(r.flags()) & 0xc0) as u8,
                    index: i as u32,
                });
            }
            keys.par_sort_unstable_by(|a, b| {
                strnum_cmp(a.name(&arena), b.name(&arena))
                    .then_with(|| a.read12.cmp(&b.read12))
                    .then_with(|| a.index.cmp(&b.index))
            });
            let order: Vec<u32> = keys.iter().map(|k| k.index).collect();
            apply_permutation(records, order);
        }
    }
}

/// Qname position in the arena, the READ1/READ2 bits, and the original index.
#[derive(Clone, Copy)]
struct NameKey {
    start: u32,
    len: u32,
    read12: u8,
    index: u32,
}

impl NameKey {
    fn name<'a>(&self, arena: &'a [u8]) -> &'a [u8] {
        let start = self.start as usize;
        &arena[start..start + self.len as usize]
    }
}

/// `(tid, pos, reverse-strand)` as htslib orders them; unplaced sorts last.
fn coord_key(rec: &RecordBuf) -> (u32, u32, u8) {
    let tid = rec
        .reference_sequence_id()
        .map(|i| i as u32)
        .unwrap_or(u32::MAX);
    let pos = rec
        .alignment_start()
        .map(|p| usize::from(p) as u32)
        .unwrap_or(u32::MAX);
    let rev = ((u16::from(rec.flags()) >> 4) & 1) as u8;
    (tid, pos, rev)
}

/// Permute `records` so position `i` holds `records[order[i]]`, rotating each
/// cycle in place. `order` is consumed: its entries double as visited marks.
fn apply_permutation(records: &mut [RecordBuf], mut order: Vec<u32>) {
    /// Marks a rotated slot; indices come from `records.len()` and never reach it.
    const DONE: u32 = u32::MAX;
    debug_assert_eq!(records.len(), order.len());
    debug_assert!(records.len() < DONE as usize);

    for start in 0..records.len() {
        if order[start] == DONE {
            continue;
        }
        let mut cur = start;
        loop {
            let src = order[cur] as usize;
            order[cur] = DONE;
            if src == start {
                break;
            }
            records.swap(cur, src);
            cur = src;
        }
    }
}

/// Coordinate-sort and optionally mark duplicates in memory, returning the sorted records.
pub fn sort_and_markdup_in_memory(
    mut records: Vec<RecordBuf>,
    markdup: bool,
) -> Result<Vec<RecordBuf>> {
    sort_records_by_kind(&mut records, SortKind::Coordinate);
    if markdup {
        let opts = crate::markdup::MarkdupOptions::default();
        let _ = mark_duplicates_in_memory(&mut records, &opts).context("mark duplicates")?;
    }
    Ok(records)
}

pub(crate) fn build_key(rec: &RecordBuf, kind: SortKind) -> SortKey {
    match kind {
        SortKind::Coordinate => {
            let tid = rec
                .reference_sequence_id()
                .map(|i| i as u32)
                .unwrap_or(u32::MAX);
            let pos = rec
                .alignment_start()
                .map(|p| usize::from(p) as u32)
                .unwrap_or(u32::MAX);
            // FLAG_REVERSE (0x10): forward (0) sorts before reverse (1), as htslib does.
            let rev = ((u16::from(rec.flags()) >> 4) & 1) as u8;
            SortKey::Coord { tid, pos, rev }
        }
        SortKind::QueryName => {
            let qname = rec.name().map(|n| n.to_vec()).unwrap_or_default();
            let read12 = (u16::from(rec.flags()) & 0xc0) as u8;
            SortKey::Name { qname, read12 }
        }
    }
}

pub fn run(args: SortArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let kind = if args.name_sort {
        SortKind::QueryName
    } else {
        SortKind::Coordinate
    };
    // `-m` is the budget for the whole sort, not per thread as in samtools.
    let mem_budget = resolve_memory_hint(
        &args.memory,
        AUTO_MIN_MEMORY_BYTES,
        AUTO_MEMORY_NUM,
        AUTO_MEMORY_DEN,
    )?;
    let chunk_capacity_hint = (mem_budget / ESTIMATED_BYTES_PER_RECORD).clamp(1024, 4 << 20);

    // mmap the input when it's on a regular file (tmp SAM from the
    // aligner is). The kernel turns the parser's reads into trivial page
    // faults; on cold files we still benefit from larger readahead.
    let open_opts = OpenOptions { mmap: true };
    let mut reader =
        BamReader::open_with_options(&args.input, args.reference.as_deref(), open_opts)
            .context("open input")?;
    let mut header = reader.header().clone();
    set_header_sort_order(&mut header, kind);
    append_pg(&mut header, &PgInfo::new("sort", !args.no_pg))?;

    let mut current_chunk: Vec<RecordBuf> = Vec::with_capacity(chunk_capacity_hint);
    let mut chunk_bytes: usize = 0;
    let mut spilled: Vec<PathBuf> = Vec::new();
    let tmpdir = match &args.tmpdir {
        Some(p) => {
            std::fs::create_dir_all(p).ok();
            TempDir::new_in(p)?
        }
        None => TempDir::new()?,
    };
    let mut spill_index = 0usize;
    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        let flags = u16::from(rec.flags());
        if args.require_flags != 0 && (flags & args.require_flags) != args.require_flags {
            continue;
        }
        if args.filter_flags != 0 && (flags & args.filter_flags) != 0 {
            continue;
        }
        chunk_bytes += approx_heap_size(&rec);
        // The reader reallocates the emptied buffers either way; cloning would
        // only add a copy on top.
        current_chunk.push(std::mem::take(&mut rec));
        if chunk_bytes >= mem_budget {
            spilled.push(spill_chunk(
                &mut current_chunk,
                &header,
                tmpdir.path(),
                spill_index,
                kind,
            )?);
            spill_index += 1;
            chunk_bytes = 0;
        }
    }

    if spilled.is_empty() {
        // Sort by key. We pre-extract keys once so par_sort doesn't allocate
        // a fresh qname Vec for every comparison (rayon's sort comparator
        // is called O(N log N) times).
        sort_records_by_kind(&mut current_chunk, kind);

        // ── Fused markdup ────────────────────────────────────────────────
        //
        // When `--markdup` is on AND the whole input fit in a single
        // in-memory chunk, we run the dup-decision pass right here on the
        // sorted slice. The standalone `kira-bam markdup` does the same
        // logic in two BAM read passes (~60s × 2 on chr20 30× WGS); this
        // path saves both passes — the records never leave RAM.
        if args.markdup {
            let opts = MarkdupOptions {
                barcode_tag: args.markdup_barcode_tag.clone(),
                ancient: args.markdup_mode_ancient,
                supplementary: args.markdup_supplementary,
            };
            let _stats =
                mark_duplicates_in_memory(&mut current_chunk, &opts).context("mark duplicates")?;
            append_pg(&mut header, &PgInfo::new("markdup", !args.no_pg))?;
        }
        let fmt = pick_fmt(&args);
        // Multi-threaded BGZF compression for BAM output. Cap at 8 workers
        // — empirically BGZF block compression saturates around there on
        // consumer SSDs and going wider just steals CPU from the rest of
        // the pipeline. Reserve 1 thread for the producer (the main one
        // issuing write_record_buf calls).
        let write_opts = WriteOptions {
            compression_workers: args.threads.saturating_sub(1).clamp(1, 8),
            compression_level: args
                .compression_level
                .and_then(bgzf::io::writer::CompressionLevel::new),
        };
        let mut writer = BamWriter::create_with_options(
            args.output.as_deref(),
            header,
            fmt,
            args.reference.as_deref(),
            write_opts,
        )?;
        writer.write_header()?;
        if matches!(fmt, OutputFormat::Bam | OutputFormat::UncompressedBam) {
            // Encode records to raw BAM bytes in parallel chunks, then stream through BGZF.
            let n_chunks = rayon::current_num_threads().max(1);
            let chunk_len = current_chunk.len().div_ceil(n_chunks).max(1);
            let encoded: Vec<(Vec<u8>, u64)> = {
                let hdr = writer.header();
                current_chunk
                    .par_chunks(chunk_len)
                    .map(|recs| crate::io::encode_records_into(hdr, recs))
                    .collect()
            };
            let skipped: u64 = encoded.iter().map(|(_, n)| n).sum();
            for (buf, _) in &encoded {
                writer.write_preencoded(buf)?;
            }
            report_skipped(skipped);
        } else {
            let mut skipped: u64 = 0;
            for r in current_chunk.iter() {
                if let Err(e) = writer.write_record_buf(r) {
                    if !is_record_level_error(&e) {
                        return Err(e);
                    }
                    skipped += 1;
                    if skipped <= 5 {
                        eprintln!("[kira-bam sort] skip record: {e:#}");
                    }
                }
            }
            report_skipped(skipped);
        }
        writer.finish()?;
        return Ok(());
    }
    if args.markdup {
        eprintln!(
            "[kira-bam sort] warning: --markdup requires single-chunk in-memory sort \
             (current memory hint produces {} spill chunks). Falling back to sort-only; \
             run `kira-bam markdup` on the output, or bump `-m`/use `-m auto`.",
            spilled.len() + (current_chunk.is_empty() as usize ^ 1)
        );
    }

    if !current_chunk.is_empty() {
        spilled.push(spill_chunk(
            &mut current_chunk,
            &header,
            tmpdir.path(),
            spill_index,
            kind,
        )?);
    }

    merge_spills(spilled, header, &args, kind)?;
    Ok(())
}

fn report_skipped(skipped: u64) {
    if skipped > 0 {
        eprintln!("[kira-bam sort] {skipped} records skipped as invalid");
    }
}

fn pick_fmt(args: &SortArgs) -> OutputFormat {
    if let Some(p) = &args.output
        && let Some(ext) = p.extension().and_then(|e| e.to_str())
    {
        if ext.eq_ignore_ascii_case("cram") {
            return OutputFormat::Cram;
        }
        if ext.eq_ignore_ascii_case("sam") {
            return OutputFormat::Sam;
        }
    }
    if args.uncompressed {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    }
}

fn set_header_sort_order(header: &mut Header, kind: SortKind) {
    use noodles_sam::header::record::value::map::header::sort_order;
    let so: &[u8] = match kind {
        SortKind::Coordinate => sort_order::COORDINATE,
        SortKind::QueryName => sort_order::QUERY_NAME,
    };
    if let Some(hdr) = header.header_mut().as_mut() {
        hdr.other_fields_mut()
            .insert(SORT_ORDER, BString::from(so.to_vec()));
    } else {
        let mut m = Map::<noodles_sam::header::record::value::map::Header>::new(Version::new(1, 6));
        m.other_fields_mut()
            .insert(SORT_ORDER, BString::from(so.to_vec()));
        *header.header_mut() = Some(m);
    }
}

fn spill_chunk(
    chunk: &mut Vec<RecordBuf>,
    header: &Header,
    dir: &std::path::Path,
    idx: usize,
    kind: SortKind,
) -> Result<PathBuf> {
    // Stable within-chunk; the k-way merge breaks equal-key ties by spill index,
    // so the external sort preserves input order as samtools does.
    sort_records_by_kind(chunk, kind);
    let path = dir.join(format!("kira-bam-spill-{idx:06}.bam"));
    let file = File::create(&path).context("create spill file")?;
    let buf = BufWriter::with_capacity(1 << 20, file);
    // FAST (level 1) instead of NONE — works around a level-0 / stored-block
    // round-trip regression in noodles-bgzf 0.47 + libdeflate.
    let bgz = bgzf::io::writer::Builder::default()
        .set_compression_level(bgzf::io::writer::CompressionLevel::FAST)
        .build_from_writer(buf);
    let mut writer = bam::io::Writer::from(bgz);
    writer.write_header(header).context("write spill header")?;
    for r in chunk.drain(..) {
        // Broken records would fail the final emit anyway; drop now to keep the run alive.
        let _ = writer.write_alignment_record(header, &r);
    }
    AlignmentWrite::finish(&mut writer, header).context("finish spill")?;
    Ok(path)
}

struct HeapItem {
    key: SortKey,
    record: RecordBuf,
    src: usize,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.src == other.src
    }
}
impl Eq for HeapItem {}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed → BinaryHeap acts as min-heap on key.
        other
            .key
            .cmp(&self.key)
            .then_with(|| other.src.cmp(&self.src))
    }
}

fn merge_spills(
    spilled: Vec<PathBuf>,
    header: Header,
    args: &SortArgs,
    kind: SortKind,
) -> Result<()> {
    let mut readers: Vec<BamReader> = spilled
        .iter()
        .map(|p| BamReader::open(p).context("reopen spill"))
        .collect::<Result<_>>()?;

    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::with_capacity(readers.len());
    for (idx, r) in readers.iter_mut().enumerate() {
        let mut rec = RecordBuf::default();
        if r.read_record_buf(&mut rec)? {
            let key = build_key(&rec, kind);
            heap.push(HeapItem {
                key,
                record: rec,
                src: idx,
            });
        }
    }

    let fmt = pick_fmt(args);
    // Same multi-threaded BGZF compression the single-chunk path uses.
    let write_opts = WriteOptions {
        compression_workers: args.threads.saturating_sub(1).clamp(1, 8),
        compression_level: args
            .compression_level
            .and_then(bgzf::io::writer::CompressionLevel::new),
    };
    let mut writer = BamWriter::create_with_options(
        args.output.as_deref(),
        header,
        fmt,
        args.reference.as_deref(),
        write_opts,
    )?;
    writer.write_header()?;

    let mut skipped: u64 = 0;
    while let Some(item) = heap.pop() {
        if let Err(e) = writer.write_record_buf(&item.record) {
            if !is_record_level_error(&e) {
                return Err(e);
            }
            skipped += 1;
            if skipped <= 5 {
                eprintln!("[kira-bam sort] skip record: {e:#}");
            }
        }
        let src = item.src;
        let mut rec = item.record;
        if readers[src].read_record_buf(&mut rec)? {
            let key = build_key(&rec, kind);
            heap.push(HeapItem {
                key,
                record: rec,
                src,
            });
        }
    }
    report_skipped(skipped);
    writer.finish()?;
    Ok(())
}

#[cfg(test)]
mod sort_order_tests {
    use super::*;
    use noodles_sam::alignment::record::Flags;

    fn rec(name: &str, tid: Option<usize>, pos: usize, flags: u16) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.name_mut() = Some(name.into());
        *r.reference_sequence_id_mut() = tid;
        *r.alignment_start_mut() = noodles_core::Position::new(pos);
        *r.flags_mut() = Flags::from(flags);
        r
    }

    fn names(records: &[RecordBuf]) -> Vec<String> {
        records
            .iter()
            .map(|r| String::from_utf8_lossy(r.name().unwrap().as_ref()).into_owned())
            .collect()
    }

    /// Reference the key-and-permute path must reproduce exactly.
    fn naive(records: &[RecordBuf], kind: SortKind) -> Vec<String> {
        let mut keyed: Vec<(SortKey, usize)> = records
            .iter()
            .enumerate()
            .map(|(i, r)| (build_key(r, kind), i))
            .collect();
        keyed.sort_by(|a, b| a.0.cmp(&b.0));
        keyed
            .iter()
            .map(|(_, i)| String::from_utf8_lossy(records[*i].name().unwrap().as_ref()).into_owned())
            .collect()
    }

    fn corpus() -> Vec<RecordBuf> {
        vec![
            rec("r10", Some(1), 100, 0),
            rec("r2", Some(0), 500, 16),
            rec("r2", Some(0), 500, 0), // same coord, forward sorts before reverse
            rec("unplaced", None, 1, 4),
            rec("r1", Some(0), 500, 0), // ties with the two above on coordinate
            rec("r100", Some(0), 1, 16),
            rec("r3", Some(1), 100, 16),
        ]
    }

    #[test]
    fn coordinate_order_matches_the_reference() {
        let mut records = corpus();
        let want = naive(&records, SortKind::Coordinate);
        sort_records_by_kind(&mut records, SortKind::Coordinate);
        assert_eq!(names(&records), want);
    }

    #[test]
    fn name_order_matches_the_reference() {
        let mut records = corpus();
        let want = naive(&records, SortKind::QueryName);
        sort_records_by_kind(&mut records, SortKind::QueryName);
        assert_eq!(names(&records), want);
    }

    /// Equal keys must keep input order.
    #[test]
    fn equal_coordinate_keys_keep_input_order() {
        let mut records = corpus();
        sort_records_by_kind(&mut records, SortKind::Coordinate);
        let got = names(&records);
        let tied: Vec<&String> = got.iter().filter(|n| *n == "r2" || *n == "r1").collect();
        // Input order among the forward-strand ties at pos 500 was r2 then r1.
        assert_eq!(tied, vec!["r2", "r1", "r2"]);
    }

    #[test]
    fn permutation_reorders_and_preserves_every_record() {
        for n in [0usize, 1, 2, 3, 5, 17, 64] {
            let records: Vec<RecordBuf> = (0..n)
                .map(|i| rec(&format!("r{i:03}"), Some(0), i + 1, 0))
                .collect();
            // Reverse order is a permutation with both 2-cycles and a fixed point.
            let order: Vec<u32> = (0..n as u32).rev().collect();
            let want: Vec<String> = names(&records).into_iter().rev().collect();
            let mut got = records;
            apply_permutation(&mut got, order);
            assert_eq!(names(&got), want, "n={n}");
        }
    }

    #[test]
    fn permutation_handles_long_cycles() {
        let n = 1000usize;
        let records: Vec<RecordBuf> = (0..n)
            .map(|i| rec(&format!("r{i:04}"), Some(0), i + 1, 0))
            .collect();
        // A single n-cycle: position i takes the record from (i + 7) % n.
        let order: Vec<u32> = (0..n).map(|i| ((i + 7) % n) as u32).collect();
        let want: Vec<String> = order
            .iter()
            .map(|&j| String::from_utf8_lossy(records[j as usize].name().unwrap().as_ref()).into_owned())
            .collect();
        let mut got = records;
        apply_permutation(&mut got, order);
        assert_eq!(names(&got), want);
    }
}

#[cfg(test)]
mod strnum_tests {
    use super::strnum_cmp;
    use std::cmp::Ordering;

    fn cmp(a: &str, b: &str) -> Ordering {
        strnum_cmp(a.as_bytes(), b.as_bytes())
    }

    #[test]
    fn numbers_compare_as_numbers() {
        assert_eq!(cmp("r2", "r10"), Ordering::Less); // 2 < 10, not byte order
        assert_eq!(cmp("r10", "r2"), Ordering::Greater);
        assert_eq!(cmp("x9y", "x10y"), Ordering::Less);
        assert_eq!(cmp("read100", "read99"), Ordering::Greater);
    }

    #[test]
    fn plain_strings_byte_order() {
        assert_eq!(cmp("alpha", "bravo"), Ordering::Less);
        assert_eq!(cmp("ab", "abc"), Ordering::Less); // prefix is smaller
        assert_eq!(cmp("abc", "abc"), Ordering::Equal);
    }

    #[test]
    fn leading_zeros_value_then_offset_tiebreak() {
        // Equal numeric value: htslib breaks the tie on how many leading zeros
        // were consumed (the one that consumed more sorts greater).
        assert_eq!(cmp("00042", "42"), Ordering::Greater);
        assert_eq!(cmp("1", "01"), Ordering::Less);
        assert_eq!(cmp("042", "042"), Ordering::Equal);
    }

    #[test]
    fn mixed_segments() {
        assert_eq!(cmp("a1b2", "a1b10"), Ordering::Less);
        assert_eq!(cmp("chr1:100", "chr1:99"), Ordering::Greater);
        assert_eq!(cmp("chr2:1", "chr10:1"), Ordering::Less); // 2 < 10
    }
}
