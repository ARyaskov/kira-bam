use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::Cigar;
use rustc_hash::FxHasher;

use crate::cli::AmpliconclipArgs;
use crate::io::{
    BamReader, BamWriter, PgInfo, append_pg, install_thread_pool, report_skipped,
    write_record_tolerant,
};
use crate::types::OutputFormat;

#[derive(Debug, Clone, Copy)]
struct Iv {
    start: u64,
    end: u64,
}

pub fn run(args: AmpliconclipArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let intervals = read_bed(&args.bed)?;

    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let mut header = reader.header().clone();
    append_pg(&mut header, &PgInfo::new("ampliconclip", !args.no_pg))?;

    let fmt = if args.uncompressed {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    };
    let mut writer = BamWriter::create(Some(&args.output), header.clone(), fmt)?;
    writer.write_header()?;

    // Contig names resolved once, not rebuilt as a String per record.
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| n.to_string())
        .collect();

    let mut rec = RecordBuf::default();
    let mut skipped = 0u64;
    while reader.read_record_buf(&mut rec)? {
        if let Some(tid) = rec.reference_sequence_id() {
            let ref_name = ref_names.get(tid).map(String::as_str).unwrap_or("");
            if let Some(ivs) = intervals.get(ref_name) {
                clip_record(&mut rec, ivs, args.hard_clip);
            }
        }
        write_record_tolerant(&mut writer, &rec, "ampliconclip", &mut skipped)?;
    }
    writer.finish()?;
    report_skipped("ampliconclip", skipped);
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
    // How many *reference* bases the overlapping primers cover at each end.
    let mut left_ref: u64 = 0;
    let mut right_ref: u64 = 0;
    for iv in ivs {
        if iv.end < start || iv.start > end {
            continue;
        }
        if iv.start <= start && iv.end >= start {
            // 5' overlap
            left_ref = left_ref.max(iv.end - start + 1);
        }
        if iv.start <= end && iv.end >= end {
            // 3' overlap
            right_ref = right_ref.max(end - iv.start + 1);
        }
    }
    if left_ref == 0 && right_ref == 0 {
        return;
    }
    apply_clip(rec, left_ref, right_ref, hard);
}

/// Result of walking one end of a CIGAR far enough to cover `want` reference
/// bases.
#[derive(Debug, Default)]
struct Trim {
    /// Reference bases actually covered (less than `want` for a short read).
    ref_bases: u64,
    /// Read bases that fall inside the trimmed span, including any existing
    /// soft clip absorbed on the way.
    read_bases: usize,
    /// Existing hard-clip length already present at this end.
    hard_bases: usize,
    /// Ops consumed whole.
    ops_consumed: usize,
    /// Bases taken off the first op *not* consumed whole — the one the clip
    /// boundary fell inside. 0 when it landed on an op edge.
    boundary_taken: usize,
}

/// Clip `left_ref`/`right_ref` reference bases off the ends of an alignment.
///
/// The clip lengths arrive in reference space (that is what a primer interval
/// measures) but a CIGAR is written in read space, and the two only coincide
/// for a gapless read with no existing clips. Converting between them is the
/// whole job here: walk the ops, count reference and read bases separately,
/// and — for the 5' end — move POS forward by exactly the reference bases
/// removed. Treating the reference count as a read offset, as this used to,
/// mis-sliced any read with an indel or a pre-existing clip and left POS
/// pointing at the primer.
fn apply_clip(rec: &mut RecordBuf, left_ref: u64, right_ref: u64, hard: bool) {
    let ops: Vec<Op> = rec.cigar().as_ref().to_vec();

    let left = walk(ops.iter().copied(), left_ref);
    let right = walk(ops.iter().rev().copied(), right_ref);

    // The two ends must not meet in the middle: a read no longer than its
    // primers is consumed entirely.
    let (lo, hi) = (left.ops_consumed, ops.len() - right.ops_consumed);
    let mut kept: Vec<Op> = if lo >= hi {
        Vec::new()
    } else {
        ops[lo..hi].to_vec()
    };
    // Whatever is left of the two boundary ops. Applying the front trim first
    // and the back trim to the result handles the case where a single op spans
    // both boundaries.
    if left.boundary_taken > 0
        && let Some(first) = kept.first_mut()
    {
        *first = Op::new(
            first.kind(),
            first.len().saturating_sub(left.boundary_taken),
        );
    }
    if right.boundary_taken > 0
        && let Some(last) = kept.last_mut()
    {
        *last = Op::new(last.kind(), last.len().saturating_sub(right.boundary_taken));
    }

    let clip_kind = if hard { Kind::HardClip } else { Kind::SoftClip };
    let mut out: Vec<Op> = Vec::with_capacity(kept.len() + 4);

    // 5' end: existing hard clip, then the new clip covering every read base
    // removed (including soft clips that were already there).
    push_op(&mut out, Kind::HardClip, left.hard_bases);
    if hard {
        push_op(&mut out, Kind::HardClip, left.read_bases);
    } else {
        push_op(&mut out, clip_kind, left.read_bases);
    }
    for op in kept {
        push_op(&mut out, op.kind(), op.len());
    }
    if hard {
        push_op(&mut out, Kind::HardClip, right.read_bases);
    } else {
        push_op(&mut out, clip_kind, right.read_bases);
    }
    push_op(&mut out, Kind::HardClip, right.hard_bases);

    let aligned = out
        .iter()
        .any(|o| o.kind().consumes_reference() && o.kind().consumes_read());
    if !aligned {
        // Nothing of the read survives outside the primers. samtools drops such
        // reads out of the alignment rather than emitting a clip-only CIGAR.
        let flags = u16::from(rec.flags()) | 0x4;
        *rec.flags_mut() = flags.into();
        *rec.cigar_mut() = Cigar::from(Vec::new());
        return;
    }

    if hard {
        // Hard clipping physically removes the bases, so SEQ/QUAL shrink by the
        // read bases taken off each end — measured in read space, which is the
        // point of `walk` reporting them separately.
        let seq = rec.sequence().as_ref();
        let from = left.read_bases.min(seq.len());
        let to = seq.len().saturating_sub(right.read_bases).max(from);
        let new_seq: Vec<u8> = seq[from..to].to_vec();
        let qual = rec.quality_scores().as_ref();
        let new_qual: Vec<u8> = qual.get(from..to.min(qual.len())).unwrap_or(&[]).to_vec();
        *rec.sequence_mut() = new_seq.into();
        *rec.quality_scores_mut() = new_qual.into();
    }

    *rec.cigar_mut() = Cigar::from(out);

    // POS moves to the first reference base still aligned.
    if left.ref_bases > 0
        && let Some(p) = rec.alignment_start()
    {
        let new_start = usize::from(p) as u64 + left.ref_bases;
        if let Ok(p) = noodles_core::Position::try_from(new_start as usize) {
            *rec.alignment_start_mut() = Some(p);
        }
    }

    // MD and NM describe the alignment we just shortened, and MC the mate's
    // pre-clip CIGAR. Recomputing them is `calmd`'s job; leaving stale values
    // behind is worse than dropping them.
    for tag in [
        Tag::MISMATCHED_POSITIONS,
        Tag::EDIT_DISTANCE,
        Tag::MATE_CIGAR,
    ] {
        rec.data_mut().remove(&tag);
    }
}

/// Walk ops from one end until `want` reference bases have been covered.
fn walk<I: Iterator<Item = Op>>(ops: I, want: u64) -> Trim {
    let mut t = Trim::default();
    for op in ops {
        let kind = op.kind();
        let len = op.len();
        if kind == Kind::HardClip {
            // Already-removed bases: they bound nothing, but the op has to be
            // carried through to the rebuilt CIGAR.
            t.hard_bases += len;
            t.ops_consumed += 1;
            continue;
        }
        if want == 0 {
            // No clipping asked for at this end, so an existing soft clip stays
            // exactly where it is.
            break;
        }
        if kind.consumes_reference() {
            if t.ref_bases >= want {
                break;
            }
            let take = ((want - t.ref_bases) as usize).min(len);
            t.ref_bases += take as u64;
            if kind.consumes_read() {
                t.read_bases += take;
            }
            if take < len {
                // The boundary falls inside this op; the remainder stays.
                t.boundary_taken = take;
                return t;
            }
            t.ops_consumed += 1;
        } else {
            // Read-only (I/S) or P. These move no reference boundary, so they
            // are absorbed into the clip — including one sitting immediately
            // against it, since a leading insertion has nothing to anchor to.
            if kind.consumes_read() {
                t.read_bases += len;
            }
            t.ops_consumed += 1;
        }
    }
    t
}

/// Append `len` bases of `kind`, merging with the previous op and dropping
/// empties, so the CIGAR never carries zero-length or adjacent duplicate ops.
fn push_op(out: &mut Vec<Op>, kind: Kind, len: usize) {
    if len == 0 {
        return;
    }
    match out.last_mut() {
        Some(prev) if prev.kind() == kind => *prev = Op::new(kind, prev.len() + len),
        _ => out.push(Op::new(kind, len)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

    fn ops(spec: &[(Kind, usize)]) -> Cigar {
        Cigar::from(spec.iter().map(|&(k, n)| Op::new(k, n)).collect::<Vec<_>>())
    }

    fn show(cigar: &Cigar) -> String {
        cigar
            .as_ref()
            .iter()
            .map(|o| {
                let c = match o.kind() {
                    Kind::Match => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::Skip => 'N',
                    Kind::SoftClip => 'S',
                    Kind::HardClip => 'H',
                    Kind::Pad => 'P',
                    Kind::SequenceMatch => '=',
                    Kind::SequenceMismatch => 'X',
                };
                format!("{}{c}", o.len())
            })
            .collect()
    }

    fn record(pos: usize, cigar: &[(Kind, usize)], read_len: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name("r")
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(pos).unwrap())
            .set_cigar(ops(cigar))
            .set_sequence(Sequence::from(vec![b'A'; read_len]))
            .set_quality_scores(QualityScores::from(vec![30u8; read_len]))
            .build()
    }

    /// A CIGAR is only valid if its read-consuming ops sum to SEQ's length.
    /// The old clip rewrote one without the other.
    fn assert_cigar_matches_seq(rec: &RecordBuf) {
        let from_cigar: usize = rec
            .cigar()
            .as_ref()
            .iter()
            .filter(|o| o.kind().consumes_read())
            .map(|o| o.len())
            .sum();
        assert_eq!(
            from_cigar,
            rec.sequence().as_ref().len(),
            "CIGAR {} does not sum to SEQ length {}",
            show(rec.cigar()),
            rec.sequence().as_ref().len()
        );
    }

    fn pos_of(rec: &RecordBuf) -> usize {
        rec.alignment_start().map(usize::from).unwrap_or(0)
    }

    #[test]
    fn soft_clips_a_gapless_read_and_moves_pos() {
        // Primer covers 100-104, i.e. the first 5 reference bases.
        let mut rec = record(100, &[(Kind::Match, 10)], 10);
        apply_clip(&mut rec, 5, 0, false);
        assert_eq!(show(rec.cigar()), "5S5M");
        assert_eq!(
            pos_of(&rec),
            105,
            "POS must skip the clipped reference bases"
        );
        assert_cigar_matches_seq(&rec);
    }

    #[test]
    fn absorbs_an_existing_soft_clip() {
        // 3S already at the 5' end: those read bases are inside the new clip,
        // so the result is one 8S — not a second clip op next to the first.
        let mut rec = record(100, &[(Kind::SoftClip, 3), (Kind::Match, 10)], 13);
        apply_clip(&mut rec, 5, 0, false);
        assert_eq!(show(rec.cigar()), "8S5M");
        assert_eq!(pos_of(&rec), 105);
        assert_cigar_matches_seq(&rec);
    }

    #[test]
    fn a_deletion_consumes_reference_but_no_read_bases() {
        // 5M2D5M spans 100..=111. Clipping 7 reference bases eats 5M and 2D,
        // which is only 5 read bases — the count the CIGAR needs.
        let mut rec = record(
            100,
            &[(Kind::Match, 5), (Kind::Deletion, 2), (Kind::Match, 5)],
            10,
        );
        apply_clip(&mut rec, 7, 0, false);
        assert_eq!(show(rec.cigar()), "5S5M");
        assert_eq!(pos_of(&rec), 107);
        assert_cigar_matches_seq(&rec);
    }

    #[test]
    fn an_insertion_at_the_boundary_is_absorbed() {
        // 5M3I5M: clipping 5 reference bases takes the 5M plus the insertion
        // that would otherwise be left dangling with nothing to anchor it.
        let mut rec = record(
            100,
            &[(Kind::Match, 5), (Kind::Insertion, 3), (Kind::Match, 5)],
            13,
        );
        apply_clip(&mut rec, 5, 0, false);
        assert_eq!(show(rec.cigar()), "8S5M");
        assert_eq!(pos_of(&rec), 105);
        assert_cigar_matches_seq(&rec);
    }

    #[test]
    fn clips_the_three_prime_end_without_moving_pos() {
        let mut rec = record(100, &[(Kind::Match, 10)], 10);
        apply_clip(&mut rec, 0, 4, false);
        assert_eq!(show(rec.cigar()), "6M4S");
        assert_eq!(pos_of(&rec), 100);
        assert_cigar_matches_seq(&rec);
    }

    #[test]
    fn clips_both_ends_of_one_op() {
        let mut rec = record(100, &[(Kind::Match, 10)], 10);
        apply_clip(&mut rec, 3, 4, false);
        assert_eq!(show(rec.cigar()), "3S3M4S");
        assert_eq!(pos_of(&rec), 103);
        assert_cigar_matches_seq(&rec);
    }

    #[test]
    fn hard_clip_trims_sequence_and_qualities() {
        let mut rec = record(100, &[(Kind::SoftClip, 3), (Kind::Match, 10)], 13);
        apply_clip(&mut rec, 5, 2, true);
        assert_eq!(show(rec.cigar()), "8H3M2H");
        assert_eq!(pos_of(&rec), 105);
        assert_eq!(rec.sequence().as_ref().len(), 3);
        assert_eq!(rec.quality_scores().as_ref().len(), 3);
        assert_cigar_matches_seq(&rec);
    }

    #[test]
    fn preserves_an_existing_hard_clip() {
        let mut rec = record(100, &[(Kind::HardClip, 4), (Kind::Match, 10)], 10);
        apply_clip(&mut rec, 5, 0, false);
        assert_eq!(show(rec.cigar()), "4H5S5M");
        assert_cigar_matches_seq(&rec);
    }

    #[test]
    fn a_fully_clipped_read_becomes_unmapped() {
        let mut rec = record(100, &[(Kind::Match, 10)], 10);
        apply_clip(&mut rec, 6, 6, false);
        assert!(
            u16::from(rec.flags()) & 0x4 != 0,
            "should be flagged unmapped"
        );
        assert_eq!(show(rec.cigar()), "");
    }

    #[test]
    fn stale_md_nm_and_mc_are_dropped() {
        use noodles_sam::alignment::record_buf::data::field::Value;
        let mut rec = record(100, &[(Kind::Match, 10)], 10);
        rec.data_mut()
            .insert(Tag::MISMATCHED_POSITIONS, Value::String("10".into()));
        rec.data_mut().insert(Tag::EDIT_DISTANCE, Value::Int32(0));
        apply_clip(&mut rec, 5, 0, false);
        assert!(rec.data().get(&Tag::MISMATCHED_POSITIONS).is_none());
        assert!(rec.data().get(&Tag::EDIT_DISTANCE).is_none());
    }
}
