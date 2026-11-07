use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use rustc_hash::FxHasher;

use crate::cli::DepthArgs;
use crate::io::{BamReader, install_thread_pool};

#[derive(Debug, Clone)]
struct DepthIv {
    start: u64,
    end: u64,
}

/// One contig's BED intervals plus a cursor. The sweep visits positions in
/// increasing order, so membership only needs the cursor advanced; testing every
/// interval per position costs `positions × intervals`.
#[derive(Debug, Default)]
pub(crate) struct BedCursor {
    ivs: Vec<DepthIv>,
    at: usize,
}

impl BedCursor {
    /// Sort and merge so the intervals are disjoint and ascending.
    fn new(mut ivs: Vec<DepthIv>) -> Self {
        ivs.sort_unstable_by_key(|iv| (iv.start, iv.end));
        let mut merged: Vec<DepthIv> = Vec::with_capacity(ivs.len());
        for iv in ivs {
            match merged.last_mut() {
                Some(last) if iv.start <= last.end + 1 => last.end = last.end.max(iv.end),
                _ => merged.push(iv),
            }
        }
        Self { ivs: merged, at: 0 }
    }

    /// Is `pos` covered? Positions must be non-decreasing between calls.
    pub(crate) fn contains(&mut self, pos: u64) -> bool {
        while self.at < self.ivs.len() && self.ivs[self.at].end < pos {
            self.at += 1;
        }
        matches!(self.ivs.get(self.at), Some(iv) if pos >= iv.start)
    }

    /// Lowest and highest position this contig's intervals touch. Lets the
    /// `-a` sweep skip straight to the restricted span instead of walking the
    /// whole reference to discard it a position at a time.
    pub(crate) fn span(&self) -> Option<(u64, u64)> {
        match (self.ivs.first(), self.ivs.last()) {
            (Some(f), Some(l)) => Some((f.start, l.end)),
            _ => None,
        }
    }
}

pub(crate) type CursorMap = HashMap<String, BedCursor, std::hash::BuildHasherDefault<FxHasher>>;

/// Parse `-r` region specs into the same cursor structure `-b` uses: both are
/// just per-contig position filters. `chr`, `chr:start` and `chr:start-end`
/// are all accepted, 1-based inclusive as samtools writes them.
pub(crate) fn parse_regions(specs: &[String], refs: &[(String, usize)]) -> Result<CursorMap> {
    let mut raw: HashMap<String, Vec<DepthIv>, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    for spec in specs {
        let (name, start, end) = match spec.rsplit_once(':') {
            // A bare contig name may itself contain ':', so only treat the tail
            // as a range when it parses as one.
            Some((name, range)) => match parse_range(range) {
                Some((s, e)) => (name.to_string(), s, e),
                None => (spec.clone(), None, None),
            },
            None => (spec.clone(), None, None),
        };
        let len = refs
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, l)| *l as u64)
            .with_context(|| format!("region {spec}: no @SQ named {name} in the header"))?;
        let start = start.unwrap_or(1).max(1);
        let end = end.unwrap_or(len).min(len);
        if start > end {
            anyhow::bail!("region {spec}: start {start} is past end {end}");
        }
        raw.entry(name).or_default().push(DepthIv { start, end });
    }
    Ok(raw
        .into_iter()
        .map(|(name, ivs)| (name, BedCursor::new(ivs)))
        .collect())
}

/// `start`, `start-end` or `start-`, all 1-based inclusive.
fn parse_range(s: &str) -> Option<(Option<u64>, Option<u64>)> {
    let s = s.replace(',', "");
    match s.split_once('-') {
        Some((a, b)) => {
            let start = if a.is_empty() {
                None
            } else {
                Some(a.parse().ok()?)
            };
            let end = if b.is_empty() {
                None
            } else {
                Some(b.parse().ok()?)
            };
            if start.is_none() && end.is_none() {
                return None;
            }
            Some((start, end))
        }
        None => Some((Some(s.parse().ok()?), None)),
    }
}

fn load_bed(path: &std::path::Path) -> Result<CursorMap> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut raw: HashMap<String, Vec<DepthIv>, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let start: u64 = parts[1].parse().unwrap_or(0);
        let end: u64 = parts[2].parse().unwrap_or(0);
        raw.entry(parts[0].to_string()).or_default().push(DepthIv {
            start: start + 1,
            end,
        });
    }
    Ok(raw
        .into_iter()
        .map(|(name, ivs)| (name, BedCursor::new(ivs)))
        .collect())
}

/// A position is reported only if it passes every restriction in force. `-b`
/// and `-r` are independent filters, so both must admit it.
fn pos_allowed(filters: &mut Filters<'_>, pos: u64) -> bool {
    filters.bed.as_deref_mut().is_none_or(|c| c.contains(pos))
        && filters
            .region
            .as_deref_mut()
            .is_none_or(|c| c.contains(pos))
}

/// The `-b` and `-r` cursors for the contig currently being swept.
#[derive(Default)]
struct Filters<'a> {
    bed: Option<&'a mut BedCursor>,
    region: Option<&'a mut BedCursor>,
}

impl Filters<'_> {
    /// The narrowest `[start, end]` the restrictions allow on this contig, or
    /// `None` when nothing is restricted.
    fn span(&self) -> Option<(u64, u64)> {
        match (
            self.bed.as_deref().and_then(|c| c.span()),
            self.region.as_deref().and_then(|c| c.span()),
        ) {
            (Some(b), Some(r)) => Some((b.0.max(r.0), b.1.min(r.1))),
            (Some(s), None) | (None, Some(s)) => Some(s),
            (None, None) => None,
        }
    }
}

pub fn run(args: DepthArgs) -> Result<()> {
    install_thread_pool(args.threads);

    // samtools spells the two levels `-a` and `-aa`/`--aa`: the first fills in
    // the zero-coverage positions of contigs that have reads, the second also
    // emits contigs that have none.
    let all_positions = args.all_positions >= 1 || args.aa;
    let all_refs = args.all_positions >= 2 || args.aa;

    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let header = reader.header().clone();
    // Resolved once rather than per record.
    let refs: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .map(|(name, sq)| (name.to_string(), usize::from(sq.length())))
        .collect();

    let mut bed = match &args.bed {
        Some(p) => Some(load_bed(p)?),
        None => None,
    };
    // `-r` used to be parsed and then ignored, silently reporting the whole
    // file. It restricts positions exactly as `-b` does, so it reuses the
    // same cursors.
    let mut regions = if args.regions.is_empty() {
        None
    } else {
        Some(parse_regions(&args.regions, &refs)?)
    };

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::with_capacity(
            1 << 20,
            File::create(p).context("create output")?,
        )),
        None => Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())),
    };

    let mut state = SweepState::default();
    let mut rec = RecordBuf::default();
    // Position of the last record accepted, to verify the input really is
    // coordinate-sorted — the sweep assumes it and produces nonsense otherwise.
    let mut last: Option<(usize, usize)> = None;
    // Header index of the first contig not yet written, for `-aa`.
    let mut next_ref = 0usize;
    // Does the current contig survive `-b`/`-r` at all?
    let mut cur_allowed = true;

    while reader.read_record_buf(&mut rec)? {
        let flags = u16::from(rec.flags());
        // samtools excludes UNMAP|SECONDARY|QCFAIL|DUP by default; counting
        // those inflated every reported depth.
        if flags & args.excl_flags != 0 {
            continue;
        }
        if args.incl_flags != 0 && (flags & args.incl_flags) != args.incl_flags {
            continue;
        }
        if let Some(mapq) = rec.mapping_quality()
            && u8::from(mapq) < args.min_mapq
        {
            continue;
        }
        let tid = match rec.reference_sequence_id() {
            Some(t) => t,
            None => continue,
        };
        let start = match rec.alignment_start() {
            Some(p) => usize::from(p),
            None => continue,
        };

        if let Some((last_tid, last_pos)) = last
            && (tid < last_tid || (tid == last_tid && start < last_pos))
        {
            let (last_name, _) = ref_entry(&refs, last_tid);
            let (name, _) = ref_entry(&refs, tid);
            anyhow::bail!(
                "input is not coordinate-sorted: {name}:{start} follows {last_name}:{last_pos}. \
                 depth sweeps positions in order, so an unordered input yields wrong depths \
                 rather than an error — sort it first (`kira-bam sort`)."
            );
        }
        last = Some((tid, start));

        // Reference positions this read actually *covers*: contiguous runs of
        // M/=/X (optionally split where base quality < -Q). Deletions (D) and
        // reference skips (N, e.g. RNA-seq introns) consume reference but place
        // no base, so they are NOT covered — counting them, as the old
        // [start,end] span did, over-reported depth across indels and inflated
        // intronic coverage to read-through levels.
        let intervals = covered_intervals(&rec, start, args.min_baseq);
        if intervals.is_empty() {
            continue;
        }

        if state.current_tid != Some(tid) {
            if let Some(prev) = state.current_tid {
                finish_contig(
                    prev,
                    &refs,
                    &mut bed,
                    &mut regions,
                    all_positions,
                    &mut state,
                    &mut out,
                )?;
            }
            if all_refs {
                // Contigs between the last one and this one carry no reads;
                // `-aa` still lists their positions at depth 0.
                emit_empty_refs(next_ref..tid, &refs, &mut bed, &mut regions, &mut out)?;
            }
            next_ref = tid + 1;

            let (name, _) = ref_entry(&refs, tid);
            cur_allowed = contig_allowed(name, &bed, &regions);
            let first_pos = if all_positions {
                filters_for(name, &mut bed, &mut regions)
                    .span()
                    .map(|(s, _)| s as usize)
                    .unwrap_or(1)
            } else {
                start
            };
            state = SweepState {
                current_tid: Some(tid),
                // Without -a we skip the leading run of zero-coverage positions;
                // with it we still start no earlier than the restricted span.
                current_pos: first_pos,
                ..Default::default()
            };
        }
        if !cur_allowed {
            continue;
        }

        let (ref_name, _) = ref_entry(&refs, tid);
        let mut filters = filters_for(ref_name, &mut bed, &mut regions);
        state.emit_until(start, ref_name, all_positions, &mut filters, &mut out)?;
        state.add_intervals(&intervals);
    }

    if let Some(prev) = state.current_tid {
        finish_contig(
            prev,
            &refs,
            &mut bed,
            &mut regions,
            all_positions,
            &mut state,
            &mut out,
        )?;
    }
    if all_refs {
        emit_empty_refs(
            next_ref..refs.len(),
            &refs,
            &mut bed,
            &mut regions,
            &mut out,
        )?;
    }
    out.flush()?;
    Ok(())
}

/// Flush the tail of `tid`'s sweep, if the restrictions admit that contig.
fn finish_contig(
    tid: usize,
    refs: &[(String, usize)],
    bed: &mut Option<CursorMap>,
    regions: &mut Option<CursorMap>,
    all_positions: bool,
    state: &mut SweepState,
    out: &mut Box<dyn Write>,
) -> Result<()> {
    let (name, len) = ref_entry(refs, tid);
    if !contig_allowed(name, bed, regions) {
        return Ok(());
    }
    let mut filters = filters_for(name, bed, regions);
    state.drain(name, len, all_positions, &mut filters, out)
}

/// Emit depth-0 rows for contigs carrying no reads at all (`-aa`).
fn emit_empty_refs(
    range: std::ops::Range<usize>,
    refs: &[(String, usize)],
    bed: &mut Option<CursorMap>,
    regions: &mut Option<CursorMap>,
    out: &mut Box<dyn Write>,
) -> Result<()> {
    for i in range {
        let Some((name, len)) = refs.get(i) else {
            continue;
        };
        if !contig_allowed(name, bed, regions) {
            continue;
        }
        let mut filters = filters_for(name, bed, regions);
        let (from, to) = filters.span().unwrap_or((1, *len as u64));
        for pos in from..=to.min(*len as u64) {
            if pos_allowed(&mut filters, pos) {
                writeln!(out, "{name}\t{pos}\t0")?;
            }
        }
    }
    Ok(())
}

/// Do the active restrictions mention this contig? A `-b`/`-r` that never names
/// it excludes it entirely — the cursor simply being absent must not read as
/// "unrestricted".
fn contig_allowed(name: &str, bed: &Option<CursorMap>, regions: &Option<CursorMap>) -> bool {
    bed.as_ref().is_none_or(|m| m.contains_key(name))
        && regions.as_ref().is_none_or(|m| m.contains_key(name))
}

fn filters_for<'a>(
    name: &str,
    bed: &'a mut Option<CursorMap>,
    regions: &'a mut Option<CursorMap>,
) -> Filters<'a> {
    Filters {
        bed: bed.as_mut().and_then(|m| m.get_mut(name)),
        region: regions.as_mut().and_then(|m| m.get_mut(name)),
    }
}

fn ref_entry(refs: &[(String, usize)], tid: usize) -> (&str, usize) {
    match refs.get(tid) {
        Some((name, len)) => (name.as_str(), *len),
        None => ("", 0),
    }
}

/// Maximal runs of reference positions actually covered by aligned read bases.
///
/// Walks the CIGAR tracking the reference position (1-based) and the read index.
/// A position is covered when the op consumes both reference and read (`M`/`=`/`X`)
/// and the corresponding base quality passes `min_baseq`. `D`/`N` break a run
/// (reference advances, no base); `I`/`S` advance only the read and are
/// transparent to reference contiguity. Returns inclusive `(start, end)` runs.
fn covered_intervals(rec: &RecordBuf, start: usize, min_baseq: u8) -> Vec<(usize, usize)> {
    let quals = rec.quality_scores().as_ref();
    let mut ref_pos = start;
    let mut read_idx = 0usize;
    let mut intervals: Vec<(usize, usize)> = Vec::new();
    let mut cur: Option<(usize, usize)> = None;

    for op in rec.cigar().as_ref().iter() {
        let kind = op.kind();
        let len = op.len();
        let cons_ref = kind.consumes_reference();
        let cons_read = kind.consumes_read();
        if cons_ref && cons_read {
            for _ in 0..len {
                let q = quals.get(read_idx).copied().unwrap_or(0);
                if min_baseq == 0 || q >= min_baseq {
                    match cur {
                        Some((a, b)) if b + 1 == ref_pos => cur = Some((a, ref_pos)),
                        Some(iv) => {
                            intervals.push(iv);
                            cur = Some((ref_pos, ref_pos));
                        }
                        None => cur = Some((ref_pos, ref_pos)),
                    }
                } else if let Some(iv) = cur.take() {
                    intervals.push(iv);
                }
                ref_pos += 1;
                read_idx += 1;
            }
        } else if cons_ref {
            // D / N: reference advances, base absent → close the current run.
            if let Some(iv) = cur.take() {
                intervals.push(iv);
            }
            ref_pos += len;
        } else if cons_read {
            // I / S: read advances, reference position unchanged (transparent).
            read_idx += len;
        }
        // H / P: consume neither; nothing to do.
    }
    if let Some(iv) = cur {
        intervals.push(iv);
    }
    intervals
}

/// Streaming per-position depth via a difference-array sweep over a min-heap of
/// `(position, delta)` events. Each covered interval `[a, b]` contributes `+1`
/// at `a` and `-1` at `b + 1`; the running sum is the depth at each position.
#[derive(Default)]
struct SweepState {
    current_tid: Option<usize>,
    /// Next reference position to emit (1-based). 0 = nothing emitted yet.
    current_pos: usize,
    cur_cov: i64,
    events: BinaryHeap<Reverse<(usize, i32)>>,
    /// Largest covered position seen, for draining without `-a`.
    max_end: usize,
}

impl SweepState {
    fn add_intervals(&mut self, intervals: &[(usize, usize)]) {
        for &(a, b) in intervals {
            self.events.push(Reverse((a, 1)));
            self.events.push(Reverse((b + 1, -1)));
            if b > self.max_end {
                self.max_end = b;
            }
        }
    }

    fn apply_events_up_to(&mut self, pos: usize) {
        while let Some(&Reverse((epos, _))) = self.events.peek() {
            if epos <= pos {
                let Reverse((_, delta)) = self.events.pop().unwrap();
                self.cur_cov += delta as i64;
            } else {
                break;
            }
        }
    }

    fn emit_until(
        &mut self,
        target: usize,
        ref_name: &str,
        all_positions: bool,
        filters: &mut Filters<'_>,
        out: &mut Box<dyn Write>,
    ) -> Result<()> {
        if self.current_pos == 0 {
            self.current_pos = 1; // positions are 1-based
        }
        while self.current_pos < target {
            let p = self.current_pos;
            self.apply_events_up_to(p);
            let cov = self.cur_cov.max(0) as usize;
            if (cov > 0 || all_positions) && pos_allowed(filters, p as u64) {
                writeln!(out, "{ref_name}\t{p}\t{cov}")?;
            }
            self.current_pos += 1;
        }
        Ok(())
    }

    fn drain(
        &mut self,
        ref_name: &str,
        ref_len: usize,
        all_positions: bool,
        filters: &mut Filters<'_>,
        out: &mut Box<dyn Write>,
    ) -> Result<()> {
        let target = if all_positions {
            // A restricted sweep stops at the end of the restriction rather
            // than walking the rest of the reference to discard it.
            let end = filters
                .span()
                .map(|(_, e)| (e as usize).min(ref_len))
                .unwrap_or(ref_len);
            end + 1
        } else {
            self.max_end + 1
        };
        self.emit_until(target, ref_name, all_positions, filters, out)
    }
}
