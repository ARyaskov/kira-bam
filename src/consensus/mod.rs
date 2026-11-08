use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};

use crate::cli::ConsensusArgs;
use crate::io::BamReader;
use crate::pileup::PileupIter;

pub fn run(args: ConsensusArgs) -> Result<()> {
    let reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let header = reader.header().clone();

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut pileup = PileupIter::new(reader, args.min_mapq, args.min_baseq);

    let mut current_tid: Option<usize> = None;
    let mut last_pos: usize = 0;
    let mut line_pos: usize = 0;

    let mut col: Vec<(u8, u8)> = Vec::new();
    while let Some((tid, pos)) = pileup.next_column_into(&mut col)? {
        if current_tid != Some(tid) {
            if current_tid.is_some() {
                writeln!(out)?;
            }
            let name = header
                .reference_sequences()
                .get_index(tid)
                .map(|(n, _)| n.to_string())
                .unwrap_or_default();
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
        if col.len() > args.max_depth as usize {
            col.truncate(args.max_depth as usize);
        }
        let base = majority(&col, args.call_fraction);
        write_base(&mut out, base, &mut line_pos, args.line_width)?;
        last_pos = pos;
    }
    if current_tid.is_some() {
        writeln!(out)?;
    }
    Ok(())
}

fn write_base(out: &mut Box<dyn Write>, b: u8, line_pos: &mut usize, width: usize) -> Result<()> {
    out.write_all(&[b])?;
    *line_pos += 1;
    if width > 0 && *line_pos >= width {
        out.write_all(b"\n")?;
        *line_pos = 0;
    }
    Ok(())
}

fn majority(col: &[(u8, u8)], threshold: f32) -> u8 {
    let mut counts = [0u32; 6]; // A C G T N other
    for &(b, _) in col {
        let i = match b.to_ascii_uppercase() {
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
    fn consensus_bases(max_depth: u32, records: &[RecordBuf]) -> String {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("in.bam");
        let fasta = dir.path().join("out.fa");
        {
            let mut w = BamWriter::create(Some(&bam), header(), OutputFormat::Bam).unwrap();
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
            max_depth,
            call_fraction: 0.5,
            line_width: 60,
            reference: None,
        };
        run(args).unwrap();

        // Skip the `>chr1` header line; collect the wrapped base calls.
        let text = std::fs::read_to_string(&fasta).unwrap();
        text.lines().skip(1).collect::<String>()
    }

    #[test]
    fn max_depth_caps_the_column_before_the_call() {
        // The first read deep calls `C`. Without the cap the column is 2 deep
        // C vs T, a tie, and this base call flipped; with `-d 1` the cap makes
        // the call follow the surviving read.
        let bases = consensus_bases(1, &[read("r1", b"CCCCC"), read("r2", b"TTTTT")]);
        assert_eq!(
            bases, "CCCCC",
            "call must reflect -d 1, not the tie: {bases:?}"
        );
    }

    #[test]
    fn zero_depth_calls_every_position_n() {
        let bases = consensus_bases(0, &[read("r1", b"ACGTA"), read("r2", b"ACGTA")]);
        assert_eq!(
            bases, "NNNNN",
            "an empty capped column must call N: {bases:?}"
        );
    }
}
