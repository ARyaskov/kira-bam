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

fn load_bed(
    path: &std::path::Path,
) -> Result<HashMap<String, Vec<DepthIv>, std::hash::BuildHasherDefault<FxHasher>>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out: HashMap<String, Vec<DepthIv>, std::hash::BuildHasherDefault<FxHasher>> =
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
        out.entry(parts[0].to_string()).or_default().push(DepthIv {
            start: start + 1,
            end,
        });
    }
    Ok(out)
}

fn pos_in_bed(ivs: Option<&Vec<DepthIv>>, pos: u64) -> bool {
    match ivs {
        Some(ivs) => ivs.iter().any(|iv| pos >= iv.start && pos <= iv.end),
        None => true,
    }
}

pub fn run(args: DepthArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let bed = match &args.bed {
        Some(p) => Some(load_bed(p)?),
        None => None,
    };

    let mut reader = BamReader::open(&args.input).context("open input")?;
    let header = reader.header().clone();

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::with_capacity(
            1 << 20,
            File::create(p).context("create output")?,
        )),
        None => Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())),
    };

    let mut state = SweepState::default();
    let mut rec = RecordBuf::default();

    while reader.read_record_buf(&mut rec)? {
        let flags = u16::from(rec.flags());
        if (flags & 0x4) != 0 {
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
                let (ref_name, ref_len) = ref_name_len(&header, prev);
                let bed_for_ref = bed.as_ref().and_then(|b| b.get(&ref_name));
                state.drain(&ref_name, ref_len, args.all_positions, bed_for_ref, &mut out)?;
            }
            state = SweepState {
                current_tid: Some(tid),
                // Without -a we skip the leading run of zero-coverage positions.
                current_pos: if args.all_positions { 0 } else { start },
                ..Default::default()
            };
        }

        let (ref_name, _) = ref_name_len(&header, tid);
        let bed_for_ref = bed.as_ref().and_then(|b| b.get(&ref_name));
        state.emit_until(start, &ref_name, args.all_positions, bed_for_ref, &mut out)?;
        state.add_intervals(&intervals);
    }

    if let Some(prev) = state.current_tid {
        let (ref_name, ref_len) = ref_name_len(&header, prev);
        let bed_for_ref = bed.as_ref().and_then(|b| b.get(&ref_name));
        state.drain(&ref_name, ref_len, args.all_positions, bed_for_ref, &mut out)?;
    }
    out.flush()?;
    Ok(())
}

fn ref_name_len(header: &noodles_sam::Header, tid: usize) -> (String, usize) {
    let name = header
        .reference_sequences()
        .get_index(tid)
        .map(|(n, _)| n.to_string())
        .unwrap_or_default();
    let len = header
        .reference_sequences()
        .get_index(tid)
        .map(|(_, sq)| usize::from(sq.length()))
        .unwrap_or(0);
    (name, len)
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
        bed: Option<&Vec<DepthIv>>,
        out: &mut Box<dyn Write>,
    ) -> Result<()> {
        if self.current_pos == 0 {
            self.current_pos = 1; // positions are 1-based
        }
        while self.current_pos < target {
            let p = self.current_pos;
            self.apply_events_up_to(p);
            let cov = self.cur_cov.max(0) as usize;
            if pos_in_bed(bed, p as u64) && (cov > 0 || all_positions) {
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
        bed: Option<&Vec<DepthIv>>,
        out: &mut Box<dyn Write>,
    ) -> Result<()> {
        let target = if all_positions {
            ref_len + 1
        } else {
            self.max_end + 1
        };
        self.emit_until(target, ref_name, all_positions, bed, out)
    }
}
