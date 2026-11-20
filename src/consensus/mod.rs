use std::io::Write;

use anyhow::{Context, Result};

use crate::cli::ConsensusArgs;
use crate::io::BamReader;
use crate::io::BufferedOutput;
use crate::pileup::{ObsKind, PileupIter, PileupObs};

pub fn run(args: ConsensusArgs) -> Result<()> {
    let reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let header = reader.header().clone();

    let mut out = BufferedOutput::create(args.output.as_deref())?;

    let mut pileup = PileupIter::new(reader, args.min_mapq, args.min_baseq);

    // @SQ order and lengths drive the output: the FASTA must list every
    // reference in header order, padded to its LN, not just the covered ones.
    let ref_seqs: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .map(|(name, sq)| (name.to_string(), usize::from(sq.length())))
        .collect();

    let mut current_tid: Option<usize> = None;
    // First contig not (partially) written yet; contigs before that either
    // had columns or were emitted whole when the sweep jumped past them.
    let mut undone: usize = 0;
    let mut last_pos: usize = 0;
    let mut line_pos: usize = 0;

    let mut col: Vec<PileupObs> = Vec::new();
    while let Some((tid, pos)) = pileup.next_column_into(&mut col)? {
        if current_tid != Some(tid) {
            if let Some(prev) = current_tid {
                // The covered span of the previous contig is complete: pad
                // from the last pileup column to its @SQ LN, numbering the
                // whole tail with N.
                while last_pos < ref_seqs[prev].1 {
                    write_base(&mut out, b'N', &mut line_pos, args.line_width)?;
                    last_pos += 1;
                }
                // `write_base` already wrapped at the width, so a full final
                // line must not gain a blank line.
                if line_pos > 0 {
                    writeln!(out)?;
                }
            }
            // Contigs the pileup skipped entirely (no records) can never
            // open a line of their own; emit them whole, in header order.
            for (name, len) in &ref_seqs[undone..tid] {
                writeln!(out, ">{name}")?;
                let mut lp = 0;
                for _ in 0..*len {
                    write_base(&mut out, b'N', &mut lp, args.line_width)?;
                }
                if lp > 0 {
                    writeln!(out)?;
                }
            }
            undone = tid + 1;
            let name = &ref_seqs[tid].0;
            writeln!(out, ">{name}")?;
            current_tid = Some(tid);
            last_pos = 0;
            line_pos = 0;
        }
        // Fill gaps with N.
        while last_pos + 1 < pos {
            write_base(&mut out, b'N', &mut line_pos, args.line_width)?;
            last_pos += 1;
        }
        // `-d` is the samtools consensus MINIMUM depth: below it the
        // position stays `N`; 0 disables the minimum entirely. Like
        // bam_consensus.c, the gate counts the column after `-Q` filtering
        // (the engine already dropped those entries), with reference skips
        // excluded (`ref_skip` entries are skipped before `depth++`) but
        // deletions still counting as observations.
        let depth = col.iter().filter(|o| o.kind != ObsKind::RefSkip).count() as u32;
        let base = if depth < args.min_depth {
            b'N'
        } else {
            majority(&col, args.call_fraction)
        };
        write_base(&mut out, base, &mut line_pos, args.line_width)?;
        last_pos = pos;
    }
    if let Some(prev) = current_tid {
        while last_pos < ref_seqs[prev].1 {
            write_base(&mut out, b'N', &mut line_pos, args.line_width)?;
            last_pos += 1;
        }
        // `write_base` already wrapped at the width, so a full final line
        // must not gain a blank line.
        if line_pos > 0 {
            writeln!(out)?;
        }
    }
    for (name, len) in &ref_seqs[undone..] {
        writeln!(out, ">{name}")?;
        let mut lp = 0;
        for _ in 0..*len {
            write_base(&mut out, b'N', &mut lp, args.line_width)?;
        }
        if lp > 0 {
            writeln!(out)?;
        }
    }
    out.finish()
}

fn write_base(out: &mut dyn Write, b: u8, line_pos: &mut usize, width: usize) -> Result<()> {
    out.write_all(&[b])?;
    *line_pos += 1;
    if width > 0 && *line_pos >= width {
        out.write_all(b"\n")?;
        *line_pos = 0;
    }
    Ok(())
}

fn majority(col: &[PileupObs], threshold: f32) -> u8 {
    let mut counts = [0u32; 6]; // A C G T N other
    for o in col {
        // Deletions (`*`) and reference skips (`>`/`<`) are depth, not
        // callable bases: only aligned bases take part in the majority.
        if o.kind != ObsKind::Base {
            continue;
        }
        let i = match o.base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 4,
            _ => 5,
        };
        counts[i] += 1;
    }
    let total: u32 = counts.iter().sum();
    if total == 0 {
        return b'N';
    }
    let (best_i, best) = counts
        .iter()
        .enumerate()
        .max_by_key(|&(_, &c)| c)
        .unwrap_or((4, &0));
    if (*best as f32) / (total as f32) < threshold {
        return b'N';
    }
    match best_i {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::BamWriter;
    use crate::types::OutputFormat;
    use noodles_core::Position;
    use noodles_sam::Header;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles_sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    fn header() -> Header {
        let len = NonZeroUsize::try_from(1_000).unwrap();
        Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(len))
            .build()
    }

    fn read(name: &str, seq: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            // Start at 1: later starts shift every call right behind N padding
            // and the assertion would have to model the padding, not the cap.
            .set_alignment_start(Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5)]))
            .set_sequence(Sequence::from(seq.to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 5]))
            .build()
    }

    /// Run consensus over a BAM of stacked reads and return the called bases.
    fn consensus_bases(min_depth: u32, records: &[RecordBuf]) -> String {
        consensus_for(header(), records, min_depth)
            .lines()
            .skip(1)
            .collect::<String>()
    }

    /// The full FASTA text: every `>` sequence line plus the wrapped calls.
    fn consensus_for(header: Header, records: &[RecordBuf], min_depth: u32) -> String {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("in.bam");
        let fasta = dir.path().join("out.fa");
        {
            let mut w = BamWriter::create(Some(&bam), header, OutputFormat::Bam).unwrap();
            w.write_header().unwrap();
            for r in records {
                w.write_record_buf(r).unwrap();
            }
            w.finish().unwrap();
        }

        let args = ConsensusArgs {
            input: bam,
            output: Some(fasta.clone()),
            min_mapq: 0,
            min_baseq: 0,
            min_depth,
            call_fraction: 0.5,
            line_width: 60,
            reference: None,
        };
        run(args).unwrap();

        std::fs::read_to_string(&fasta).unwrap()
    }

    #[test]
    fn depth_below_min_depth_stays_n() {
        // Two stacked reads: `-d 3` is above the observed depth 2, so every
        // covered position stays N; `-d 2` calls the majority.
        let reads = [read("r1", b"CCCCC"), read("r2", b"CCCCC")];
        assert_eq!(
            consensus_bases(3, &reads),
            format!("{}{}", "N".repeat(5), "N".repeat(995)),
            "observed depth 2 below -d 3 must stay N",
        );
        assert_eq!(
            consensus_bases(2, &reads),
            format!("{}{}", "C".repeat(5), "N".repeat(995)),
            "observed depth 2 meets -d 2: {reads:?}",
        );
    }

    #[test]
    fn zero_min_depth_never_gates_the_call() {
        // `-d 0` disables the minimum: single coverage still calls its base.
        let bases = consensus_bases(0, &[read("r1", b"ACGTA")]);
        assert_eq!(
            bases,
            format!("{}{}", "ACGTA", "N".repeat(995)),
            "with no minimum a 1-deep column must call: {bases:?}",
        );
    }

    #[test]
    fn deletions_and_skips_are_not_callable_bases() {
        // r1 covers 1..=5 with C; r2 deletes reference base 3 and skips
        // reference base 4. The column depths still count them (2 and 2), so
        // `-d 1` passes the minimum — but the majority must come from the
        // aligned base alone, not from `*`/`>` observations.
        let r2 = RecordBuf::builder()
            .set_name("r2")
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Skip, 1),
                Op::new(Kind::Match, 1),
            ]))
            .set_sequence(Sequence::from(b"CCC".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 3]))
            .build();
        let bases = consensus_bases(1, &[read("r1", b"CCCCC"), r2]);
        assert_eq!(
            bases,
            format!("{}{}", "CCCCC", "N".repeat(995)),
            "deletion and skip observations must not outvote the aligned base",
        );
    }

    #[test]
    fn refskip_excluded_from_min_depth_but_deletions_count() {
        // r1 covers 1..=5 with C at depth 1; r2 deletes reference base 3 and
        // skips reference base 4. With `-d 2`: position 3 counts the
        // deletion (depth 2, call C), position 4 loses the refskip
        // (depth 1 < 2, stays N) — bam_consensus.c skips ref_skip entries
        // before `depth++` but deletions still count.
        let r2 = RecordBuf::builder()
            .set_name("r2")
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Skip, 1),
                Op::new(Kind::Match, 1),
            ]))
            .set_sequence(Sequence::from(b"CCT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 3]))
            .build();
        let bases = consensus_bases(2, &[read("r1", b"CCCCC"), r2]);
        // The discriminating calls: 3 = C (deletion counts, gate passes) and
        // 4 = N (refskip excluded, depth 1 < 2). Position 5 is a C/T tie at
        // fraction 0.5, which `--call-fraction 0.5` calls on the later
        // maximum (T).
        assert_eq!(
            bases,
            format!("CCCNT{}", "N".repeat(995)),
            "deletions count toward `-d`, refskips do not: {bases:?}",
        );
    }

    #[test]
    fn output_extends_to_the_reference_length() {
        // chr1 is 8 long; a 5-base read covers 1..=5, so 6..=8 must be N.
        let len = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8).unwrap()),
            )
            .build();
        let text = consensus_for(len, &[read("r1", b"CCCCC")], 1);
        assert_eq!(text, ">chr1\nCCCCCNNN\n");
    }

    #[test]
    fn sequences_filling_the_width_exactly_leave_no_blank_line() {
        // Uncovered LN 60 (== line_width): one full wrapped all-N line and
        // nothing after it.
        let all_n = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(60).unwrap()),
            )
            .build();
        let text = consensus_for(all_n, &[], 1);
        assert_eq!(text, format!(">chr1\n{}\n", "N".repeat(60)));

        // Covered with LN 65 (multiple of width plus 5 calls): the wrapped
        // call+tail line ends exactly at width 60 and 5 tail Ns follow on
        // their own line, with no blank line between.
        let covered = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(65).unwrap()),
            )
            .build();
        let text = consensus_for(covered, &[read("r1", b"CCCCC")], 1);
        let covered_line = format!("{}{}", "C".repeat(5), "N".repeat(55));
        assert_eq!(text, format!(">chr1\n{covered_line}\n{}\n", "N".repeat(5)));
    }

    #[test]
    fn uncovered_contig_is_emitted_as_all_n() {
        // chr2 has no records: it must still appear, in header order, all N.
        let header = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(10).unwrap()),
            )
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(6).unwrap()),
            )
            .build();
        let text = consensus_for(header, &[read("r1", b"CCCCC")], 1);
        assert_eq!(text, ">chr1\nCCCCCNNNNN\n>chr2\nNNNNNN\n");
    }

    #[test]
    fn empty_input_emits_every_reference() {
        // A small contig with no records at all is emitted whole and filled
        // to its @SQ LN with N (larger contigs only change the wrapping).
        let len = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8).unwrap()),
            )
            .build();
        let text = consensus_for(len, &[], 1);
        assert_eq!(text, ">chr1\nNNNNNNNN\n");
    }
}
