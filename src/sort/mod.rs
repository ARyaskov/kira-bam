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
    resolve_memory_hint,
};
use crate::markdup::{MarkdupOptions, mark_duplicates_in_memory};
use crate::types::OutputFormat;

// Conservative per-record heap estimate used to turn a memory budget into a
// chunk size. A `RecordBuf` owns its name, sequence, quality scores, CIGAR and
// the tag data map, so a typical 150 bp paired read costs ~400–600 B live, not
// the 256 B assumed previously — which let `-m auto` (3/4 of RAM) overshoot peak
// RSS by 2–4×. 512 keeps the in-memory chunk safely within the stated budget;
// worst case we spill a little sooner.
const ESTIMATED_BYTES_PER_RECORD: usize = 512;

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

/// Sort the records in `records` by the requested key, in parallel.
///
/// The naïve `par_sort_unstable_by(|a, b| build_key(a, kind).cmp(&build_key(b, kind)))`
/// allocates a qname `Vec<u8>` twice per comparison — for 12 M records sorted
/// by coordinate that's ~480 M Vec allocations and dominates the sort itself.
///
/// Pattern: extract `(SortKey, RecordBuf)` pairs once in parallel, sort them
/// by key (which doesn't re-run `build_key`), then re-collect the records.
/// Each `build_key` invocation runs exactly once per record.
fn sort_records_by_kind(records: &mut Vec<RecordBuf>, kind: SortKind) {
    let original = std::mem::take(records);
    let mut keyed: Vec<(SortKey, RecordBuf)> = original
        .into_par_iter()
        .map(|r| (build_key(&r, kind), r))
        .collect();
    // STABLE sort: equal keys keep their input order, matching samtools' tie-break.
    keyed.par_sort_by(|a, b| a.0.cmp(&b.0));
    records.reserve_exact(keyed.len());
    records.extend(keyed.into_iter().map(|(_, r)| r));
}

/// Coordinate-sort and optionally mark duplicates in memory, returning the sorted records.
pub fn sort_and_markdup_in_memory(
    mut records: Vec<RecordBuf>,
    markdup: bool,
) -> Result<Vec<RecordBuf>> {
    sort_records_by_kind(&mut records, SortKind::Coordinate);
    if markdup {
        let opts = crate::markdup::MarkdupOptions {
            barcode_tag: None,
            ancient: false,
        };
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
    let mem_per_thread = resolve_memory_hint(
        &args.memory,
        AUTO_MIN_MEMORY_BYTES,
        AUTO_MEMORY_NUM,
        AUTO_MEMORY_DEN,
    )?;
    let chunk_records = (mem_per_thread / ESTIMATED_BYTES_PER_RECORD).max(1024);

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

    let mut current_chunk: Vec<RecordBuf> = Vec::with_capacity(chunk_records);
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
        current_chunk.push(rec.clone());
        if current_chunk.len() >= chunk_records {
            spilled.push(spill_chunk(
                &mut current_chunk,
                &header,
                tmpdir.path(),
                spill_index,
                kind,
            )?);
            spill_index += 1;
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
            let encoded: Vec<Vec<u8>> = {
                let hdr = writer.header();
                current_chunk
                    .par_chunks(chunk_len)
                    .map(|recs| crate::io::encode_records_into(hdr, recs))
                    .collect()
            };
            for buf in &encoded {
                writer.write_preencoded(buf)?;
            }
        } else {
            let mut skipped: u64 = 0;
            for r in current_chunk.iter() {
                if let Err(e) = writer.write_record_buf(r) {
                    skipped += 1;
                    if skipped <= 5 {
                        eprintln!("[kira-bam sort] skip record: {e:#}");
                    }
                }
            }
            if skipped > 0 {
                eprintln!("[kira-bam sort] {skipped} records skipped due to write errors");
            }
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
    // STABLE within-chunk sort; the k-way merge then breaks equal-key ties by
    // spill index, so the whole external sort preserves input order like samtools.
    chunk.par_sort_by(|a, b| build_key(a, kind).cmp(&build_key(b, kind)));
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
    let mut writer = BamWriter::create_with_reference(
        args.output.as_deref(),
        header,
        fmt,
        args.reference.as_deref(),
    )?;
    writer.write_header()?;

    let mut skipped: u64 = 0;
    while let Some(item) = heap.pop() {
        if let Err(e) = writer.write_record_buf(&item.record) {
            skipped += 1;
            if skipped <= 5 {
                eprintln!("[kira-bam sort] skip record: {e:#}");
            }
        }
        let src = item.src;
        let mut rec = RecordBuf::default();
        if readers[src].read_record_buf(&mut rec)? {
            let key = build_key(&rec, kind);
            heap.push(HeapItem {
                key,
                record: rec,
                src,
            });
        }
    }
    if skipped > 0 {
        eprintln!("[kira-bam sort] {skipped} records skipped due to write errors");
    }
    writer.finish()?;
    Ok(())
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
