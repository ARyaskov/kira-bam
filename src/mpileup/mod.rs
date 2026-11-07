use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};

use crate::cli::MpileupArgs;
use crate::depth::{CursorMap, parse_regions};
use crate::faidx::RefCache;
use crate::io::BamReader;
use crate::pileup::PileupIter;

pub fn run(args: MpileupArgs) -> Result<()> {
    if args.inputs.len() > 1 {
        eprintln!("[kira-bam mpileup] note: only the first BAM is processed in this version");
    }
    let input = &args.inputs[0];
    let reader =
        BamReader::open_with_reference(input, args.reference.as_deref()).context("open input")?;
    let header = reader.header().clone();

    // Indexed, one contig cached at a time — a whole-genome FASTA no longer
    // has to fit in RAM to call a pileup.
    let mut refs = match &args.reference {
        Some(p) => Some(RefCache::open(p)?),
        None => None,
    };
    // Contig names and lengths resolved once; the region parser validates
    // against the lengths.
    let ref_lens: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .map(|(name, sq)| (name.to_string(), usize::from(sq.length())))
        .collect();

    // `-r` restricts every row — covered or zero-filled — exactly as `depth -r`
    // does (same parser: overlaps merged, 1-based inclusive bounds). A contig
    // no region names stays silent, even under `-a`.
    let mut regions = if args.regions.is_empty() {
        None
    } else {
        Some(parse_regions(&args.regions, &ref_lens)?)
    };

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let all = args.all_positions;
    let mut pileup = PileupIter::new(reader, args.min_mapq, args.min_baseq);
    // One column buffer for the whole scan.
    let mut col: Vec<(u8, u8)> = Vec::new();
    // Reference the sweep is on, whether it can emit anything at all, the last
    // position it has accounted for, and the last admissible position (`-a`
    // fills trailing zeros up to it).
    let mut cur_tid: Option<usize> = None;
    let mut cur_open = false;
    let mut upto = 0u64;
    let mut span_end = 0u64;

    while let Some((tid, pos)) = pileup.next_column_into(&mut col)? {
        if cur_tid != Some(tid) {
            // The previous reference still owes trailing zero rows under `-a`.
            if all && cur_open {
                let name = ref_lens[cur_tid.unwrap()].0.as_str();
                zero_fill(&mut out, &mut refs, &mut regions, name, upto + 1, span_end)?;
            }
            let Some((name, len)) = ref_lens.get(tid) else {
                continue;
            };
            let Some((start, end)) = fill_span(&regions, name, *len) else {
                // No region names this reference: it stays silent, even
                // under `-a`.
                cur_open = false;
                cur_tid = Some(tid);
                continue;
            };
            cur_open = true;
            upto = start.saturating_sub(1);
            span_end = end;
            cur_tid = Some(tid);
        }
        let name = ref_lens[tid].0.as_str();

        if all && cur_open {
            // Zero-coverage positions between the last row and this column.
            zero_fill(
                &mut out,
                &mut refs,
                &mut regions,
                name,
                upto + 1,
                (pos as u64).saturating_sub(1),
            )?;
        }
        if cur_open && region_allows(&mut regions, name, pos as u64) {
            if col.len() > args.max_depth as usize {
                col.truncate(args.max_depth as usize);
            }
            let ref_base = match &mut refs {
                Some(cache) => cache
                    .sequence(name)?
                    .and_then(|s| s.get(pos.saturating_sub(1)))
                    .copied()
                    .unwrap_or(b'N'),
                None => b'N',
            };
            write!(out, "{name}\t{pos}\t{}\t{}", ref_base as char, col.len())?;
            // Read bases column with ref-base style markup (. for match, ACGT for mismatch).
            write!(out, "\t")?;
            for &(b, _q) in &col {
                let same = b.eq_ignore_ascii_case(&ref_base);
                if same {
                    out.write_all(b".")?;
                } else {
                    out.write_all(&[b.to_ascii_uppercase()])?;
                }
            }
            // Quality column.
            write!(out, "\t")?;
            for &(_b, q) in &col {
                // Saturating: the 0xFF "quality unavailable" sentinel overflows a
                // plain `+ 33`.
                out.write_all(&[q.saturating_add(33).min(126)])?;
            }
            writeln!(out)?;
        }
        upto = pos as u64;
    }

    if all && cur_open {
        let name = ref_lens[cur_tid.unwrap()].0.as_str();
        zero_fill(&mut out, &mut refs, &mut regions, name, upto + 1, span_end)?;
    }
    out.flush()?;
    Ok(())
}

/// With `-a`, the `[start, end]` span a reference can produce rows in: the
/// whole 1-based reference, or the `-r` span (already clamped and merged by
/// `parse_regions`). `None` when no region names the reference — nothing may
/// come from it.
fn fill_span(regions: &Option<CursorMap>, name: &str, len: usize) -> Option<(u64, u64)> {
    match regions {
        None => Some((1, len as u64)),
        Some(m) => m.get(name).and_then(|c| c.span()),
    }
}

/// Does `pos` pass the `-r` filter? Positions within one reference must come
/// non-decreasing (`BedCursor` sweeps its cursor forward only).
fn region_allows(regions: &mut Option<CursorMap>, name: &str, pos: u64) -> bool {
    regions
        .as_mut()
        .is_none_or(|m| m.get_mut(name).is_none_or(|c| c.contains(pos)))
}

/// Emit zero-depth mpileup rows for `from..=to`, keeping only admissible
/// positions. `*` placeholders keep the six-column mpileup shape, matching
/// samtools' zero-coverage rows. Never emits already-emitted positions.
fn zero_fill(
    out: &mut Box<dyn Write>,
    refs: &mut Option<RefCache>,
    regions: &mut Option<CursorMap>,
    name: &str,
    from: u64,
    to: u64,
) -> Result<()> {
    for pos in from..=to {
        if !region_allows(regions, name, pos) {
            continue;
        }
        let ref_base = match refs {
            Some(cache) => cache
                .sequence(name)?
                .and_then(|s| s.get(pos as usize - 1))
                .copied()
                .unwrap_or(b'N'),
            None => b'N',
        };
        writeln!(out, "{name}\t{pos}\t{}\t0\t*\t*", ref_base as char)?;
    }
    Ok(())
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
            .add_reference_sequence("chr2", Map::<ReferenceSequence>::new(len))
            .build()
    }

    /// A 5 bp perfect match on `(tid, start)`.
    fn rec(name: &str, tid: usize, start: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(tid)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5)]))
            .set_sequence(Sequence::from(b"ACGTA".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 5]))
            .build()
    }

    fn write_bam(path: &std::path::Path, records: &[RecordBuf]) {
        let mut w = BamWriter::create(Some(path), header(), OutputFormat::Bam).unwrap();
        w.write_header().unwrap();
        for r in records {
            w.write_record_buf(r).unwrap();
        }
        w.finish().unwrap();
    }

    /// Write the test BAM once and run mpileup over it with the given `-r`/`-a`
    /// options, returning every row's tab-separated fields.
    fn rows(records: &[RecordBuf], regions: &[&str], all: bool) -> Vec<Vec<String>> {
        let dir = tempfile::tempdir().unwrap();
        write_bam(dir.path().join("in.bam").as_path(), records);
        let out = dir.path().join("out.tsv");
        run(MpileupArgs {
            inputs: vec![dir.path().join("in.bam")],
            output: Some(out.clone()),
            reference: None,
            regions: regions.iter().map(|s| s.to_string()).collect(),
            min_mapq: 0,
            min_baseq: 0,
            max_depth: 8000,
            all_positions: all,
        })
        .unwrap();
        std::fs::read_to_string(out)
            .unwrap()
            .lines()
            .map(|l| l.split('\t').map(String::from).collect())
            .collect()
    }

    #[test]
    fn region_restricts_covered_rows() {
        let rows = rows(&[rec("r1", 0, 10)], &["chr1:8-20"], false);
        let pos: Vec<u64> = rows.iter().map(|r| r[1].parse().unwrap()).collect();
        assert_eq!(pos, vec![10, 11, 12, 13, 14]);
        assert!(rows.iter().all(|r| r[0] == "chr1" && r[3] == "1"), "{rows:?}");
    }

    #[test]
    fn all_positions_fills_region_gaps_and_unions_overlaps() {
        // The union of both regions is discharged once, in order: 1-5 zero,
        // 8-9 zero, then the 5 covered bases, then 15-20 zero.
        let rows = rows(
            &[rec("r1", 0, 10)],
            &["chr1:1-5", "chr1:8-20"],
            true,
        );
        let depths: Vec<u64> = rows.iter().map(|r| r[3].parse().unwrap()).collect();
        assert_eq!(
            depths,
            vec![0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        );
        assert!(rows.iter().all(|r| r[0] == "chr1"), "{rows:?}");
    }

    #[test]
    fn all_positions_fills_represented_references_only() {
        // Both reads sit on chr1; chr2 has none, so `-a` must not emit it.
        let rows = rows(&[rec("r1", 0, 10), rec("r2", 0, 900)], &[], true);

        // Reads to 1000-long chr1: the whole contig is filled, gaps at zero.
        assert_eq!(rows.len(), 1_000);
        assert_eq!(rows[0][3], "0");
        assert_eq!(rows[9][3], "1");
        assert_eq!(rows[899][3], "1");
        assert_eq!(rows[910][3], "0");
        assert!(rows.iter().all(|r| r[0] == "chr1"), "{rows:?}");
    }
}
