use std::collections::HashMap;
use std::collections::hash_map::Entry;

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::data::field::Tag;
use rustc_hash::FxHasher;

use crate::cli::MarkdupArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

const FLAG_PAIRED: u16 = 0x1;
const FLAG_UNMAP: u16 = 0x4;
const FLAG_MUNMAP: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct ReadKey {
    tid: i32,
    pos5: i32,
    strand: u8,
    barcode: Vec<u8>,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct PairKey {
    a: ReadKey,
    b: ReadKey,
}

impl PairKey {
    fn new(a: ReadKey, b: ReadKey) -> Self {
        if (a.tid, a.pos5, a.strand, a.barcode.as_slice())
            <= (b.tid, b.pos5, b.strand, b.barcode.as_slice())
        {
            PairKey { a, b }
        } else {
            PairKey { a: b, b: a }
        }
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
enum DupKey {
    Single(ReadKey),
    Pair(PairKey),
}

#[derive(Default, Debug)]
pub struct DupStats {
    pub examined: u64,
    pub paired_examined: u64,
    pub single_examined: u64,
    pub duplicate_pair: u64,
    pub duplicate_single: u64,
    pub optical_duplicates: u64,
}

impl DupStats {
    pub fn total_dups(&self) -> u64 {
        self.duplicate_pair + self.duplicate_single
    }
}

/// The best template seen so far for a given [`DupKey`]. We track templates (by
/// qname), not individual reads: both mates of a pair hash to the same key, so a
/// per-read "winner" would mark one mate of the kept pair as a duplicate. The
/// winner's qname lets pass 2 mark *all* reads of every losing template while
/// keeping *both* reads of the winner.
#[derive(Clone, Debug)]
struct TemplateWinner {
    score: u32,
    qname: Vec<u8>,
}

/// Per-record dedup facts extracted once and reused in both passes.
struct Examined {
    key: DupKey,
    is_pair: bool,
    qname: Vec<u8>,
    score: u32,
}

pub fn run(args: MarkdupArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let optical_distance = args.optical_distance;

    // Pass 1: scan, build dedup decisions without holding records in RAM.
    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input (pass 1)")?;
    let mut header = reader.header().clone();
    append_pg(&mut header, &PgInfo::new("markdup", !args.no_pg))?;

    let barcode_tag = parse_barcode_tag(args.barcode_tag.as_deref())?;
    let mut winners: HashMap<DupKey, TemplateWinner, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    // Optical-duplicate accounting is opt-in (`-d > 0`): only then do we retain
    // per-group flow-cell coordinates, so the default path stays lean on memory.
    let mut optical: Option<HashMap<DupKey, OpticalGroup, std::hash::BuildHasherDefault<FxHasher>>> =
        (optical_distance > 0).then(HashMap::default);
    let mut stats = DupStats::default();

    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        if let Some(ex) = examine(&rec, barcode_tag, args.mode_ancient) {
            stats.examined += 1;
            if ex.is_pair {
                stats.paired_examined += 1;
            } else {
                stats.single_examined += 1;
            }
            if let Some(groups) = optical.as_mut()
                && let Some(loc) = parse_optical_loc(rec.name().map(AsRef::as_ref).unwrap_or(&[]))
            {
                optical_collect(groups, &ex, loc);
            }
            update_winner(&mut winners, ex);
        }
    }

    // Tally optical duplicates: the subset of PCR duplicates explained by
    // flow-cell proximity (Picard-style — per group, size minus #clusters).
    if let Some(groups) = &optical {
        stats.optical_duplicates = groups
            .values()
            .map(|g| count_optical(&g.locs, optical_distance) * if g.is_pair { 2 } else { 1 })
            .sum();
    }

    // Pass 2: stream, mark every read of a losing template, write.
    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input (pass 2)")?;
    let fmt = if args
        .output
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("cram"))
        .unwrap_or(false)
    {
        OutputFormat::Cram
    } else if args.uncompressed {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    };
    let mut writer = BamWriter::create_with_reference(
        Some(&args.output),
        header,
        fmt,
        args.reference.as_deref(),
    )?;
    writer.write_header()?;

    let mut rec = RecordBuf::default();
    let mut read_total: u64 = 0;
    let mut emitted: u64 = 0;
    while reader.read_record_buf(&mut rec)? {
        read_total += 1;
        let dup = match examine(&rec, barcode_tag, args.mode_ancient) {
            Some(ex) => is_template_dup(&winners, &ex, &mut stats),
            None => false,
        };
        if dup {
            let flags = u16::from(rec.flags()) | FLAG_DUP;
            *rec.flags_mut() = noodles_sam::alignment::record::Flags::from(flags);
        }
        if !(args.remove && dup) {
            writer.write_record_buf(&rec)?;
            emitted += 1;
        }
    }
    writer.finish()?;

    if let Some(stats_path) = args.stats {
        write_stats(
            &stats_path,
            &stats,
            read_total,
            emitted,
            optical_distance > 0,
        )
        .context("write stats")?;
    }
    Ok(())
}

/// Options for the in-memory markdup pass (`mark_duplicates_in_memory`).
#[derive(Clone, Debug, Default)]
pub struct MarkdupOptions {
    /// SAM tag (2 chars) holding the cell barcode for single-cell dedup;
    /// `None` means barcode-less standard dedup.
    pub barcode_tag: Option<String>,
    /// "Ancient DNA" mode: ignore strand in the dup key (matches samtools
    /// `--mode s`). Off by default; on, more aggressive dedup.
    pub ancient: bool,
}

/// Mark PCR/library duplicates on `records` in place, returning stats.
///
/// This in-memory helper (used by the fused `sort --markdup` path) does not
/// compute optical duplicates — `DupStats::optical_duplicates` stays 0. The
/// standalone `markdup -d` command does the optical pass (see [`run`]).
///
/// **Preconditions:**
/// * `records` is sorted by coordinate (TID, position).
/// * Each record's mate pointer (PNEXT, mate-reference-id) is correct, and ideally
///   the `MC` tag is present (run `kira-bam fixmate -m` first) so the mate's
///   *unclipped* 5′ end can be reconstructed for clipped reads. Without `MC` we
///   key on the clipped PNEXT, which is an approximation for soft/hard-clipped
///   mates.
///
/// **Algorithm:** the same template-aware logic as the standalone two-pass
/// markdup (`markdup::run`). We extract a `DupKey` per record, pick the
/// best-scoring *template* (by qname) per key, then flag every read of the
/// losing templates (FLAG 0x400) — keeping both mates of the winner. The
/// advantage over the two-pass version is that records never leave RAM, so we
/// save one full BAM read+decode pass (~60 s for chr20 30× WGS).
///
/// The pass is parallelisable over chunks because the dup key is
/// position-local; future work could shard by (TID, position-bucket) and
/// merge winner maps. Today we keep it single-threaded for correctness
/// margin.
pub fn mark_duplicates_in_memory(
    records: &mut [RecordBuf],
    options: &MarkdupOptions,
) -> Result<DupStats> {
    use noodles_sam::alignment::record::Flags;
    let barcode_tag = parse_barcode_tag(options.barcode_tag.as_deref())?;
    let ancient = options.ancient;

    let mut winners: HashMap<DupKey, TemplateWinner, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    let mut stats = DupStats::default();

    // Pass A — pick the winning template per key.
    for rec in records.iter() {
        if let Some(ex) = examine(rec, barcode_tag, ancient) {
            stats.examined += 1;
            if ex.is_pair {
                stats.paired_examined += 1;
            } else {
                stats.single_examined += 1;
            }
            update_winner(&mut winners, ex);
        }
    }

    // Pass B — flag every read whose template lost (both mates of a losing pair,
    // never a mate of the winner).
    for rec in records.iter_mut() {
        let dup = match examine(rec, barcode_tag, ancient) {
            Some(ex) => is_template_dup(&winners, &ex, &mut stats),
            None => false,
        };
        if dup {
            let cur = u16::from(rec.flags());
            *rec.flags_mut() = Flags::from(cur | FLAG_DUP);
        }
    }
    Ok(stats)
}

/// Extract the dedup facts for one record, or `None` if it isn't a primary
/// mapped read (and thus never examined for duplication).
fn examine(rec: &RecordBuf, barcode_tag: Option<Tag>, ancient: bool) -> Option<Examined> {
    let flags = u16::from(rec.flags());
    if (flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_UNMAP)) != 0 {
        return None;
    }
    let rk = read_key(rec, barcode_tag, ancient)?;
    let (key, is_pair) = match mate_key(rec, barcode_tag, ancient) {
        Some(mk) => (DupKey::Pair(PairKey::new(rk, mk)), true),
        None => (DupKey::Single(rk), false),
    };
    let qname: Vec<u8> = rec.name().map(|n| n.to_vec()).unwrap_or_default();
    // Template score: a pair scores own quality + the mate's quality sum (the `ms`
    // tag from `fixmate -m`), so BOTH mates compute the same value and we rank
    // whole templates. Without `ms` we use the read's own quality — still
    // template-consistent for marking, just a coarser tiebreak.
    let score =
        base_quality_sum(rec) + if is_pair { mate_score(rec).unwrap_or(0) } else { 0 };
    Some(Examined {
        key,
        is_pair,
        qname,
        score,
    })
}

/// Fold one examined record into the per-key winner map (higher score wins;
/// qname breaks ties deterministically). Seeing the winner's own mate again is a
/// no-op.
fn update_winner(
    winners: &mut HashMap<DupKey, TemplateWinner, std::hash::BuildHasherDefault<FxHasher>>,
    ex: Examined,
) {
    match winners.entry(ex.key) {
        Entry::Vacant(v) => {
            v.insert(TemplateWinner {
                score: ex.score,
                qname: ex.qname,
            });
        }
        Entry::Occupied(mut o) => {
            let w = o.get_mut();
            if ex.qname == w.qname {
                return; // same template (the other mate) — not a competitor
            }
            if ex.score > w.score || (ex.score == w.score && ex.qname < w.qname) {
                w.score = ex.score;
                w.qname = ex.qname;
            }
        }
    }
}

/// Is this record's template a duplicate (i.e. not the winner for its key)?
/// Updates the duplicate counters on the way.
fn is_template_dup(
    winners: &HashMap<DupKey, TemplateWinner, std::hash::BuildHasherDefault<FxHasher>>,
    ex: &Examined,
    stats: &mut DupStats,
) -> bool {
    match winners.get(&ex.key) {
        Some(w) if w.qname != ex.qname => {
            if ex.is_pair {
                stats.duplicate_pair += 1;
            } else {
                stats.duplicate_single += 1;
            }
            true
        }
        _ => false,
    }
}

/// Mate score (`ms` tag, written by `fixmate -m`): the sum of the mate's base
/// qualities. Used so paired templates are ranked by combined quality.
fn mate_score(rec: &RecordBuf) -> Option<u32> {
    use noodles_sam::alignment::record_buf::data::field::Value;
    match rec.data().get(&Tag::new(b'm', b's')) {
        Some(Value::Int32(v)) => Some((*v).max(0) as u32),
        Some(Value::UInt32(v)) => Some(*v),
        Some(Value::Int16(v)) => Some((*v).max(0) as u32),
        Some(Value::UInt16(v)) => Some(*v as u32),
        Some(Value::Int8(v)) => Some((*v).max(0) as u32),
        Some(Value::UInt8(v)) => Some(*v as u32),
        _ => None,
    }
}

// ── Optical-duplicate detection ───────────────────────────────────────────────
//
// Two PCR duplicates are *optical* when they come from the same flow-cell tile
// and lie within `-d` pixels of each other (they're one cluster the sequencer
// read twice, not an independent library molecule). We parse the cluster
// coordinate out of the Illumina read name, then within each duplicate group
// union-find reads that are tile-local and within distance; the optical count is
// `group_size - cluster_count` (Picard semantics).

/// Cluster coordinate parsed from an Illumina read name.
#[derive(Clone, Copy, Debug)]
struct OpticalLoc {
    /// Hash of the tile identity (everything up to `x`: instrument…lane:tile).
    tile: u64,
    x: i32,
    y: i32,
}

/// Per-`DupKey` optical bookkeeping (only populated when `-d > 0`).
struct OpticalGroup {
    is_pair: bool,
    /// qname hashes already recorded, so the two mates of one template count once.
    seen: Vec<u64>,
    locs: Vec<OpticalLoc>,
}

fn fxhash_bytes(b: &[u8]) -> u64 {
    use std::hash::Hasher;
    let mut h = FxHasher::default();
    h.write(b);
    h.finish()
}

/// Parse an Illumina CASAVA-style read name `…:tile:x:y`. Takes the last two
/// colon-separated fields as `x`/`y` and hashes the prefix as the tile identity.
/// Returns `None` if the name doesn't look like a coordinate name (too few
/// fields, or non-numeric `x`/`y`) — such reads simply aren't optical-eligible.
fn parse_optical_loc(name: &[u8]) -> Option<OpticalLoc> {
    let name = match name.iter().position(|&b| b == b' ') {
        Some(i) => &name[..i],
        None => name,
    };
    let colons: Vec<usize> = name
        .iter()
        .enumerate()
        .filter(|&(_, &b)| b == b':')
        .map(|(i, _)| i)
        .collect();
    if colons.len() < 4 {
        return None; // need >= 5 fields to resemble instrument:run:fc:lane:tile:x:y
    }
    let last = colons[colons.len() - 1];
    let second_last = colons[colons.len() - 2];
    let x = parse_u32_ascii(&name[second_last + 1..last])? as i32;
    let y = parse_u32_ascii(&name[last + 1..])? as i32;
    let tile = fxhash_bytes(&name[..second_last]);
    Some(OpticalLoc { tile, x, y })
}

fn parse_u32_ascii(b: &[u8]) -> Option<u32> {
    if b.is_empty() {
        return None;
    }
    let mut v: u32 = 0;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((c - b'0') as u32)?;
    }
    Some(v)
}

/// Record one template's cluster location under its dup key (deduping the two
/// mates of a pair by qname hash).
fn optical_collect(
    groups: &mut HashMap<DupKey, OpticalGroup, std::hash::BuildHasherDefault<FxHasher>>,
    ex: &Examined,
    loc: OpticalLoc,
) {
    let h = fxhash_bytes(&ex.qname);
    let g = groups.entry(ex.key.clone()).or_insert_with(|| OpticalGroup {
        is_pair: ex.is_pair,
        seen: Vec::new(),
        locs: Vec::new(),
    });
    if !g.seen.contains(&h) {
        g.seen.push(h);
        g.locs.push(loc);
    }
}

/// Number of optical duplicates among a duplicate group: `size - cluster_count`,
/// where two members cluster when they share a tile and are within `d` pixels in
/// both axes (union-find).
fn count_optical(locs: &[OpticalLoc], d: u32) -> u64 {
    let n = locs.len();
    if n < 2 {
        return 0;
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        let mut root = i;
        while parent[root] != root {
            root = parent[root];
        }
        let mut cur = i;
        while parent[cur] != root {
            let next = parent[cur];
            parent[cur] = root;
            cur = next;
        }
        root
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if locs[i].tile == locs[j].tile
                && (locs[i].x - locs[j].x).unsigned_abs() <= d
                && (locs[i].y - locs[j].y).unsigned_abs() <= d
            {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    let clusters = (0..n).filter(|&i| find(&mut parent, i) == i).count();
    (n - clusters) as u64
}

fn parse_barcode_tag(s: Option<&str>) -> Result<Option<Tag>> {
    let Some(s) = s else { return Ok(None) };
    if s.len() != 2 {
        anyhow::bail!("barcode tag must be 2 ASCII characters, got {s:?}");
    }
    let bytes = s.as_bytes();
    Ok(Some(Tag::from([bytes[0], bytes[1]])))
}

fn extract_barcode(rec: &RecordBuf, tag: Option<Tag>) -> Vec<u8> {
    use noodles_sam::alignment::record_buf::data::field::Value;
    let Some(t) = tag else { return Vec::new() };
    match rec.data().get(&t) {
        Some(Value::String(s)) => s.to_vec(),
        Some(Value::Hex(s)) => s.to_vec(),
        _ => Vec::new(),
    }
}

fn unclipped_5_prime(rec: &RecordBuf) -> i32 {
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record_buf::Cigar;
    let flags = u16::from(rec.flags());
    let pos = rec
        .alignment_start()
        .map(|p| usize::from(p) as i32)
        .unwrap_or(0);
    let cigar: &Cigar = rec.cigar();
    if (flags & FLAG_REVERSE) == 0 {
        // Leading-clip length on forward strand.
        let mut clip: i32 = 0;
        for op in cigar.as_ref().iter() {
            match op.kind() {
                Kind::SoftClip | Kind::HardClip => clip += op.len() as i32,
                _ => break,
            }
        }
        pos - clip
    } else {
        // Reverse strand: compute pos + ref_consumed - 1 + trailing_clip.
        //
        // Old code did `let ops: Vec<_> = cigar.as_ref().iter().collect()` and
        // then walked it twice (rev for trailing clip, forward for ref
        // consumed). That allocates a Vec per record — 6 M allocations on
        // a chr20 30× WGS run, half the records being reverse-strand. Here
        // we walk the CIGAR a single time, accumulating both quantities, and
        // do the "first non-clip from the end" check by latching: any
        // non-clip op resets `trailing_clip` to 0. The result is identical
        // because by definition only the contiguous trailing clip ops
        // contribute.
        let mut ref_consumed: i32 = 0;
        let mut trailing_clip: i32 = 0;
        for op in cigar.as_ref().iter() {
            if op.kind().consumes_reference() {
                ref_consumed += op.len() as i32;
            }
            match op.kind() {
                Kind::SoftClip | Kind::HardClip => trailing_clip += op.len() as i32,
                _ => trailing_clip = 0,
            }
        }
        pos + ref_consumed - 1 + trailing_clip
    }
}

fn read_key(rec: &RecordBuf, barcode_tag: Option<Tag>, ancient: bool) -> Option<ReadKey> {
    let flags = u16::from(rec.flags());
    if (flags & FLAG_UNMAP) != 0 {
        return None;
    }
    let tid = rec.reference_sequence_id()? as i32;
    let strand = if ancient {
        0
    } else if (flags & FLAG_REVERSE) != 0 {
        1
    } else {
        0
    };
    Some(ReadKey {
        tid,
        pos5: unclipped_5_prime(rec),
        strand,
        barcode: extract_barcode(rec, barcode_tag),
    })
}

fn mate_key(rec: &RecordBuf, barcode_tag: Option<Tag>, ancient: bool) -> Option<ReadKey> {
    let flags = u16::from(rec.flags());
    if (flags & FLAG_PAIRED) == 0 || (flags & FLAG_MUNMAP) != 0 {
        return None;
    }
    let tid = rec.mate_reference_sequence_id()? as i32;
    let pnext = rec.mate_alignment_start().map(|p| usize::from(p) as i32)?;
    let mate_reverse = (flags & FLAG_MATE_REVERSE) != 0;
    // Duplicate detection must key on the mate's *unclipped* 5′ end, the same
    // coordinate we compute for the read itself in `unclipped_5_prime`. PNEXT is
    // the mate's clipped leftmost mapped base, so soft/hard-clipped mates would
    // otherwise be mis-grouped (two PCR copies with different clipping get
    // different keys, or two distinct fragments collapse). The mate's CIGAR
    // lives in the MC tag (written by `fixmate -m`); when present we reconstruct
    // its unclipped 5′ exactly. Without MC we fall back to PNEXT and accept the
    // approximation — run `kira-bam fixmate -m` first for samtools-grade dedup.
    let pos5 = mate_unclipped_5_prime(pnext, mate_reverse, mate_cigar_bytes(rec).as_deref());
    let strand = if ancient {
        0
    } else if mate_reverse {
        1
    } else {
        0
    };
    Some(ReadKey {
        tid,
        pos5,
        strand,
        barcode: extract_barcode(rec, barcode_tag),
    })
}

/// Raw bytes of the mate CIGAR (MC tag), if present.
fn mate_cigar_bytes(rec: &RecordBuf) -> Option<Vec<u8>> {
    use noodles_sam::alignment::record_buf::data::field::Value;
    match rec.data().get(&Tag::MATE_CIGAR) {
        Some(Value::String(s)) => Some(s.to_vec()),
        _ => None,
    }
}

/// Mate's unclipped 5′ reference coordinate, derived from PNEXT, the mate's
/// strand, and its CIGAR (MC tag). Mirrors [`unclipped_5_prime`] for the read's
/// own end. Falls back to the clipped PNEXT when MC is absent.
fn mate_unclipped_5_prime(pnext: i32, mate_reverse: bool, mc: Option<&[u8]>) -> i32 {
    match mc.and_then(parse_mate_cigar_clips) {
        Some((leading_clip, trailing_clip, ref_span)) => {
            if !mate_reverse {
                pnext - leading_clip
            } else {
                pnext + ref_span - 1 + trailing_clip
            }
        }
        None => pnext,
    }
}

/// Parse a CIGAR string (e.g. `5S100M2S`) into
/// `(leading_clip, trailing_clip, reference_span)`, counting both soft and hard
/// clips at each end. Returns `None` for an empty/`*` CIGAR.
fn parse_mate_cigar_clips(mc: &[u8]) -> Option<(i32, i32, i32)> {
    let mut num: i64 = 0;
    let mut leading: i32 = 0;
    let mut trailing: i32 = 0;
    let mut ref_span: i32 = 0;
    let mut seen_nonclip = false;
    let mut saw_op = false;
    for &c in mc {
        if c.is_ascii_digit() {
            num = num * 10 + (c - b'0') as i64;
        } else {
            let len = num as i32;
            num = 0;
            saw_op = true;
            match c {
                b'S' | b'H' => {
                    if !seen_nonclip {
                        leading += len;
                    }
                    // Tentatively trailing; reset by the next non-clip op.
                    trailing += len;
                }
                b'M' | b'D' | b'N' | b'=' | b'X' => {
                    seen_nonclip = true;
                    trailing = 0;
                    ref_span += len;
                }
                _ => {
                    // I / P consume no reference and aren't clips.
                    seen_nonclip = true;
                    trailing = 0;
                }
            }
        }
    }
    if !saw_op { None } else { Some((leading, trailing, ref_span)) }
}

fn base_quality_sum(rec: &RecordBuf) -> u32 {
    base_quality_sum_bytes(rec.quality_scores().as_ref())
}

/// `Σ q[i]` over bytes `q[i] >= 15`. Picarrd's quality threshold (the
/// "quality of useful bases") that samtools markdup uses to pick the
/// best read in a dup group.
///
/// On x86_64 with AVX2, this runs at ~1 byte/cycle via:
///   1. 32-byte unsigned-compare against threshold 14 → 0xFF/0x00 mask,
///   2. AND the original bytes with the mask to zero out below-threshold lanes,
///   3. `sad_epu8` against zero to fold 32 lanes into 4 × u64 partial sums,
///   4. extract & sum.
///
/// vs the scalar baseline at ~3 byte/cycle in iterator form (filter + map +
/// sum on a Vec<u8>). For 12 M records × ~150 bytes that's ~1.5 s scalar →
/// ~0.3 s vector on a 4 GHz Zen 3 — modest but free.
#[inline]
fn base_quality_sum_bytes(q: &[u8]) -> u32 {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: runtime check.
            return unsafe { base_quality_sum_avx2(q) };
        }
    }
    base_quality_sum_scalar(q)
}

#[inline]
fn base_quality_sum_scalar(q: &[u8]) -> u32 {
    let mut total: u32 = 0;
    for &x in q {
        if x >= 15 {
            total += x as u32;
        }
    }
    total
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn base_quality_sum_avx2(q: &[u8]) -> u32 {
    use std::arch::x86_64::{
        __m256i, _mm_add_epi64, _mm_cvtsi128_si64, _mm_unpackhi_epi64, _mm256_add_epi64,
        _mm256_and_si256, _mm256_cmpgt_epi8, _mm256_extracti128_si256, _mm256_loadu_si256,
        _mm256_sad_epu8, _mm256_set1_epi8, _mm256_setzero_si256, _mm256_xor_si256,
    };

    // We must match `base_quality_sum_scalar` *exactly*, including for bytes
    // >= 128. The naive `_mm256_cmpgt_epi8(v, 14)` does a SIGNED compare, so the
    // `0xFF` missing-quality sentinel (-1 as i8) is dropped, while the scalar
    // path keeps it (`0xFF >= 15`). That divergence made markdup's winner choice
    // depend on whether the host CPU has AVX2 — a reproducibility hazard. Here we
    // emulate an UNSIGNED compare `q > 14` by biasing both operands by 0x80
    // before the signed `cmpgt`, so AVX2 and scalar agree on every byte 0..=255.
    let n = q.len();
    let mut i: usize = 0;
    let bias = _mm256_set1_epi8(0x80u8 as i8);
    // (14 ^ 0x80) is the biased threshold; `(q ^ 0x80) > (14 ^ 0x80)` ⇔ unsigned q > 14.
    let threshold = _mm256_set1_epi8((14u8 ^ 0x80) as i8);
    let zero = _mm256_setzero_si256();
    let mut acc = _mm256_setzero_si256();
    while i + 32 <= n {
        // SAFETY: loop bound.
        let v = unsafe { _mm256_loadu_si256(q.as_ptr().add(i) as *const __m256i) };
        let mask = _mm256_cmpgt_epi8(_mm256_xor_si256(v, bias), threshold);
        // AND against the *original* bytes so SAD sums true values (incl. 255).
        let masked = _mm256_and_si256(v, mask);
        let sad = _mm256_sad_epu8(masked, zero); // 4 × u64 partial sums
        acc = _mm256_add_epi64(acc, sad);
        i += 32;
    }
    // Reduce the 4 × u64 lanes.
    let lo = _mm256_extracti128_si256::<0>(acc);
    let hi = _mm256_extracti128_si256::<1>(acc);
    let sum128 = _mm_add_epi64(lo, hi);
    let lo64 = _mm_cvtsi128_si64(sum128) as u64;
    let hi64 = _mm_cvtsi128_si64(_mm_unpackhi_epi64(sum128, sum128)) as u64;
    let mut total: u32 = (lo64 + hi64) as u32;
    // Tail.
    for &x in &q[i..] {
        if x >= 15 {
            total += x as u32;
        }
    }
    total
}

fn write_stats(
    path: &std::path::Path,
    stats: &DupStats,
    read_total: u64,
    written: u64,
    optical_computed: bool,
) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "COMMAND: kira-bam markdup")?;
    writeln!(f, "READ {read_total}")?;
    writeln!(f, "WRITTEN {written}")?;
    writeln!(f, "EXCLUDED {}", read_total.saturating_sub(stats.examined))?;
    writeln!(f, "EXAMINED {}", stats.examined)?;
    writeln!(f, "PAIRED {}", stats.paired_examined)?;
    writeln!(f, "SINGLE {}", stats.single_examined)?;
    writeln!(f, "DUPLICATE PAIR {}", stats.duplicate_pair)?;
    writeln!(f, "DUPLICATE SINGLE {}", stats.duplicate_single)?;
    writeln!(f, "DUPLICATE TOTAL {}", stats.total_dups())?;
    // Optical duplicates are computed only when the user opts in with `-d`.
    // Without it we report NA rather than a fake `0` that downstream tooling
    // would read as "0 optical duplicates found".
    if optical_computed {
        writeln!(f, "DUPLICATE OPTICAL {}", stats.optical_duplicates)?;
    } else {
        writeln!(f, "DUPLICATE OPTICAL NA (not computed; pass -d to enable)")?;
    }
    Ok(())
}

#[cfg(test)]
mod base_quality_tests {
    use super::*;

    #[test]
    fn empty() {
        assert_eq!(base_quality_sum_bytes(&[]), 0);
    }

    #[test]
    fn below_threshold_excluded() {
        let q = vec![0, 5, 10, 14];
        assert_eq!(base_quality_sum_bytes(&q), 0);
    }

    #[test]
    fn threshold_inclusive() {
        let q = vec![15];
        assert_eq!(base_quality_sum_bytes(&q), 15);
    }

    #[test]
    fn matches_scalar_for_realistic_qualities() {
        // 150-byte read, typical Illumina quality range
        let q: Vec<u8> = (0..150u8).map(|i| (i % 41) + 20).collect();
        let scalar = base_quality_sum_scalar(&q);
        let actual = base_quality_sum_bytes(&q);
        assert_eq!(scalar, actual);
    }

    #[test]
    fn matches_scalar_for_long_random_quality_strings() {
        // Exercise tail handling. Length not a multiple of 32.
        let q: Vec<u8> = (0..255u32).map(|i| ((i * 37 + 1) % 60) as u8).collect();
        let scalar = base_quality_sum_scalar(&q);
        let actual = base_quality_sum_bytes(&q);
        assert_eq!(scalar, actual);
    }

    #[test]
    fn tail_only() {
        let q = vec![20u8; 10]; // < 32, tail-only path
        assert_eq!(base_quality_sum_bytes(&q), 200);
    }

    #[test]
    fn sentinel_0xff_matches_scalar() {
        // The 0xFF "no quality" sentinel is the historical AVX2/scalar divergence:
        // signed cmpgt dropped it, scalar kept it. Pin the behaviour explicitly.
        let q = vec![0xFFu8; 40]; // > 32 so it exercises the vector body + tail
        assert_eq!(base_quality_sum_scalar(&q), base_quality_sum_bytes(&q));
        assert_eq!(base_quality_sum_bytes(&q), 40 * 255);
    }

    proptest::proptest! {
        /// The vectorised path must equal the scalar reference for *every* byte
        /// value 0..=255 and length, so markdup's winner selection is identical
        /// regardless of whether the host CPU has AVX2.
        #[test]
        fn simd_matches_scalar_over_full_u8_domain(
            q in proptest::collection::vec(proptest::num::u8::ANY, 0..600usize)
        ) {
            proptest::prop_assert_eq!(
                base_quality_sum_scalar(&q),
                base_quality_sum_bytes(&q)
            );
        }
    }
}

#[cfg(test)]
mod optical_tests {
    use super::*;

    #[test]
    fn parses_casava_name() {
        let loc = parse_optical_loc(b"EAS139:136:FC706VJ:2:2104:15343:197393").unwrap();
        assert_eq!(loc.x, 15343);
        assert_eq!(loc.y, 197393);
    }

    #[test]
    fn rejects_non_coordinate_name() {
        assert!(parse_optical_loc(b"read1").is_none());
        assert!(parse_optical_loc(b"a:b:c").is_none()); // too few fields
        assert!(parse_optical_loc(b"i:r:fc:1:2104:xx:10").is_none()); // non-numeric x
    }

    #[test]
    fn same_tile_differs_from_other_tile() {
        let a = parse_optical_loc(b"I:1:FC:1:1101:100:100").unwrap();
        let b = parse_optical_loc(b"I:1:FC:1:1101:105:103").unwrap();
        let c = parse_optical_loc(b"I:1:FC:1:2222:100:100").unwrap(); // different tile
        assert_eq!(a.tile, b.tile);
        assert_ne!(a.tile, c.tile);
    }

    #[test]
    fn clusters_within_distance() {
        // a,b within 5px on the same tile → 1 cluster of 2 → 1 optical dup.
        // c is far away → its own cluster, 0 optical.
        let a = parse_optical_loc(b"I:1:FC:1:1101:100:100").unwrap();
        let b = parse_optical_loc(b"I:1:FC:1:1101:103:102").unwrap();
        let c = parse_optical_loc(b"I:1:FC:1:1101:9000:9000").unwrap();
        assert_eq!(count_optical(&[a, b, c], 5), 1); // 3 locs - 2 clusters
        assert_eq!(count_optical(&[a, b], 5), 1);
        assert_eq!(count_optical(&[a, c], 5), 0); // 2 locs - 2 clusters
        assert_eq!(count_optical(&[a], 5), 0);
    }

    #[test]
    fn distance_is_inclusive_and_per_axis() {
        let a = parse_optical_loc(b"I:1:FC:1:1101:100:100").unwrap();
        let on_edge = parse_optical_loc(b"I:1:FC:1:1101:105:100").unwrap(); // dx=5
        let just_over = parse_optical_loc(b"I:1:FC:1:1101:106:100").unwrap(); // dx=6
        assert_eq!(count_optical(&[a, on_edge], 5), 1);
        assert_eq!(count_optical(&[a, just_over], 5), 0);
    }
}
