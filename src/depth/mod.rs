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
        let end = match rec.alignment_end() {
            Some(p) => usize::from(p),
            None => continue,
        };
        if !overlapping_consumed_bases(&rec) {
            continue;
        }

        if state.current_tid != Some(tid) {
            // Flush previous reference fully.
            if let Some(prev) = state.current_tid {
                let ref_len = header
                    .reference_sequences()
                    .get_index(prev)
                    .map(|(_, sq)| usize::from(sq.length()))
                    .unwrap_or(0);
                let ref_name = header
                    .reference_sequences()
                    .get_index(prev)
                    .map(|(n, _)| n.to_string())
                    .unwrap_or_default();
                let bed_for_ref = bed.as_ref().and_then(|b| b.get(&ref_name));
                state.drain_to(
                    &ref_name,
                    ref_len,
                    args.all_positions,
                    bed_for_ref,
                    &mut out,
                )?;
            }
            state = SweepState {
                current_tid: Some(tid),
                ..Default::default()
            };
        }

        let ref_name = header
            .reference_sequences()
            .get_index(tid)
            .map(|(n, _)| n.to_string())
            .unwrap_or_default();
        let bed_for_ref = bed.as_ref().and_then(|b| b.get(&ref_name));
        state.advance_to(&ref_name, start, args.all_positions, bed_for_ref, &mut out)?;
        state.add_read(end);
    }

    if let Some(prev) = state.current_tid {
        let ref_name = header
            .reference_sequences()
            .get_index(prev)
            .map(|(n, _)| n.to_string())
            .unwrap_or_default();
        let ref_len = header
            .reference_sequences()
            .get_index(prev)
            .map(|(_, sq)| usize::from(sq.length()))
            .unwrap_or(0);
        let bed_for_ref = bed.as_ref().and_then(|b| b.get(&ref_name));
        state.drain_to(
            &ref_name,
            ref_len,
            args.all_positions,
            bed_for_ref,
            &mut out,
        )?;
    }
    out.flush()?;
    Ok(())
}

fn overlapping_consumed_bases(rec: &RecordBuf) -> bool {
    rec.cigar()
        .as_ref()
        .iter()
        .any(|op| op.kind().consumes_reference())
}

#[derive(Default)]
struct SweepState {
    current_tid: Option<usize>,
    current_pos: usize,
    ends: BinaryHeap<Reverse<usize>>,
}

impl SweepState {
    fn add_read(&mut self, end: usize) {
        self.ends.push(Reverse(end));
    }

    fn advance_to(
        &mut self,
        ref_name: &str,
        new_pos: usize,
        all_positions: bool,
        bed: Option<&Vec<DepthIv>>,
        out: &mut Box<dyn Write>,
    ) -> Result<()> {
        if new_pos == 0 {
            return Ok(());
        }
        let from = self.current_pos.max(1);
        for pos in from..new_pos {
            // Pop ends < pos
            while let Some(&Reverse(e)) = self.ends.peek() {
                if e < pos {
                    self.ends.pop();
                } else {
                    break;
                }
            }
            let cov = self.ends.len();
            if !pos_in_bed(bed, pos as u64) {
                continue;
            }
            if cov > 0 || all_positions {
                writeln!(out, "{ref_name}\t{pos}\t{cov}")?;
            }
        }
        self.current_pos = new_pos;
        Ok(())
    }

    fn drain_to(
        &mut self,
        ref_name: &str,
        max_pos: usize,
        all_positions: bool,
        bed: Option<&Vec<DepthIv>>,
        out: &mut Box<dyn Write>,
    ) -> Result<()> {
        let limit = if all_positions {
            max_pos
        } else {
            // Drain until heap empty.
            self.ends
                .iter()
                .map(|&Reverse(e)| e)
                .max()
                .unwrap_or(self.current_pos)
        };
        self.advance_to(ref_name, limit + 1, all_positions, bed, out)
    }
}
