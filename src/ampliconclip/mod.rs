use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Result};
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record_buf::Cigar;
use noodles_sam::alignment::RecordBuf;
use rustc_hash::FxHasher;

use crate::cli::AmpliconclipArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

#[derive(Debug, Clone, Copy)]
struct Iv {
    start: u64,
    end: u64,
}

pub fn run(args: AmpliconclipArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let intervals = read_bed(&args.bed)?;

    let mut reader = BamReader::open(&args.input).context("open input")?;
    let mut header = reader.header().clone();
    append_pg(&mut header, &PgInfo::new("ampliconclip", !args.no_pg))?;

    let fmt = if args.uncompressed {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    };
    let mut writer = BamWriter::create(Some(&args.output), header.clone(), fmt)?;
    writer.write_header()?;

    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        if let Some(tid) = rec.reference_sequence_id() {
            let ref_name: String = header
                .reference_sequences()
                .get_index(tid)
                .map(|(n, _)| n.to_string())
                .unwrap_or_default();
            if let Some(ivs) = intervals.get(&ref_name) {
                clip_record(&mut rec, ivs, args.hard_clip);
            }
        }
        let _ = writer.write_record_buf(&rec);
    }
    writer.finish()?;
    Ok(())
}

fn read_bed(
    path: &std::path::Path,
) -> Result<HashMap<String, Vec<Iv>, std::hash::BuildHasherDefault<FxHasher>>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out: HashMap<String, Vec<Iv>, std::hash::BuildHasherDefault<FxHasher>> =
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
        out.entry(parts[0].to_string()).or_default().push(Iv {
            start: start + 1,
            end,
        });
    }
    Ok(out)
}

fn clip_record(rec: &mut RecordBuf, ivs: &[Iv], hard: bool) {
    let start = match rec.alignment_start() {
        Some(p) => usize::from(p) as u64,
        None => return,
    };
    let end = match rec.alignment_end() {
        Some(p) => usize::from(p) as u64,
        None => return,
    };
    // Find leftmost primer that overlaps the read's 5' side.
    let mut left_clip: u32 = 0;
    let mut right_clip: u32 = 0;
    for iv in ivs {
        if iv.end < start || iv.start > end {
            continue;
        }
        if iv.start <= start && iv.end >= start {
            // 5' overlap
            left_clip = left_clip.max((iv.end - start + 1) as u32);
        }
        if iv.start <= end && iv.end >= end {
            // 3' overlap
            right_clip = right_clip.max((end - iv.start + 1) as u32);
        }
    }
    if left_clip == 0 && right_clip == 0 {
        return;
    }
    let new_cigar = apply_clip(rec.cigar(), left_clip, right_clip, hard);
    *rec.cigar_mut() = new_cigar;
    if hard {
        // Hard-clip also trims sequence + qual.
        let seq = rec.sequence().as_ref();
        let qual = rec.quality_scores().as_ref();
        let read_len = seq.len();
        let new_start = (left_clip as usize).min(read_len);
        let new_end = read_len.saturating_sub(right_clip as usize);
        if new_start < new_end {
            let new_seq: Vec<u8> = seq[new_start..new_end].to_vec();
            let new_qual: Vec<u8> = qual.get(new_start..new_end).unwrap_or(&[]).to_vec();
            *rec.sequence_mut() = new_seq.into();
            *rec.quality_scores_mut() = new_qual.into();
        }
    }
}

fn apply_clip(cigar: &Cigar, left: u32, right: u32, hard: bool) -> Cigar {
    let kind = if hard { Kind::HardClip } else { Kind::SoftClip };
    let ops: Vec<Op> = cigar.as_ref().iter().copied().collect();
    let mut out: Vec<Op> = Vec::new();
    if left > 0 {
        out.push(Op::new(kind, left as usize));
    }
    let mut left_remaining = left;
    let mut right_remaining = right;
    let total_len: u32 = ops.iter().map(|o| o.len() as u32).sum();
    let mut consumed: u32 = 0;
    for op in ops {
        let mut len = op.len() as u32;
        // Skip from left.
        if left_remaining > 0 && op.kind().consumes_read() {
            let take = left_remaining.min(len);
            len -= take;
            left_remaining -= take;
            if len == 0 {
                consumed += op.len() as u32;
                continue;
            }
        }
        // Determine if we'd cross right boundary.
        let pos_after = consumed + op.len() as u32;
        let right_start = total_len.saturating_sub(right_remaining);
        if right_remaining > 0 && pos_after > right_start {
            let keep = right_start.saturating_sub(consumed);
            if keep > 0 {
                out.push(Op::new(op.kind(), keep as usize));
            }
            right_remaining = 0;
        } else {
            out.push(Op::new(op.kind(), len as usize));
        }
        consumed += op.len() as u32;
    }
    if right > 0 {
        out.push(Op::new(kind, right as usize));
    }
    Cigar::from(out)
}
