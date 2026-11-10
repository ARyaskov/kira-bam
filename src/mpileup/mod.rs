use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::cli::MpileupArgs;
use crate::depth::{CursorMap, parse_regions};
use crate::faidx::RefCache;
use crate::io::{BamReader, RecordSource, indexed::IndexedInput};
use crate::merge::{RgRenames, build_tid_maps, is_sorted_compatible, merge_headers};
use crate::pileup::PileupIter;
use crate::sort::SortKind;
use noodles_sam::Header;

pub fn run(args: MpileupArgs) -> Result<()> {
    let inputs = &args.inputs;
    let readers = inputs
        .iter()
        .map(|p| {
            BamReader::open_with_reference(p, args.reference.as_deref())
                .with_context(|| format!("open {}", p.display()))
        })
        .collect::<Result<Vec<_>>>()?;
    let input_headers: Vec<Header> = readers.iter().map(|r| r.header().clone()).collect();

    // One shared @SQ dictionary = the union across inputs, in first-seen order
    // (length conflicts error), so every column is keyed by a common tid. The
    // pileup only needs that dictionary, not the read-group/program merging the
    // helper also does.
    let (union_header, _): (Header, Vec<RgRenames>) = merge_headers(
        input_headers.iter().collect::<Vec<_>>(),
        SortKind::Coordinate,
    )?;
    // Per input: local tid → union tid, matched by @SQ name.
    let tid_maps = build_tid_maps(&input_headers, &union_header);
    for (i, map) in tid_maps.iter().enumerate() {
        if map.iter().any(|m| m.is_none()) {
            bail!(
                "input {}: a @SQ name is missing from the merged dictionary",
                i
            );
        }
        if !is_sorted_compatible(map) {
            bail!(
                "input {} ({}) has @SQ in an order incompatible with the merged \
                 dictionary; the reference-column sweep cannot reorder it. \
                 Reheader it first.",
                i,
                inputs[i].display()
            );
        }
    }

    // Indexed, one contig cached at a time — a whole-genome FASTA no longer
    // has to fit in RAM to call a pileup.
    let mut refs = match &args.reference {
        Some(p) => Some(RefCache::open(p)?),
        None => None,
    };
    // Contig names and lengths resolved once; the region parser validates
    // against the lengths.
    let ref_lens: Vec<(String, usize)> = union_header
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
    let n = readers.len();
    // Eligible inputs are read through index random access: only the contigs
    // some region names — and every record of those contigs — reach the
    // sweep. Everything else (no index, SAM/CRAM, or no regions at all) keeps
    // the streaming reader. A present-but-unparsable index aborts here
    // instead of quietly scanning the file with unknown rows.
    let mut iters: Vec<PileupIter> = Vec::with_capacity(n);
    for (i, r) in readers.into_iter().enumerate() {
        let source: Box<dyn RecordSource> = match &regions {
            Some(regions) => match IndexedInput::open(&inputs[i], regions)? {
                Some(idx) => Box::new(idx),
                None => Box::new(r),
            },
            None => Box::new(r),
        };
        iters.push(PileupIter::new(source, args.min_mapq, args.min_baseq));
    }
    // Column buffer per input, reused across positions.
    let mut cols: Vec<Vec<(u8, u8)>> = vec![Vec::new(); n];
    // The position each input currently stands on, in union tid space.
    let mut heads: Vec<Option<(usize, usize)>> = Vec::with_capacity(n);
    for i in 0..n {
        heads.push(next_col(
            &mut iters[i],
            &tid_maps[i],
            &inputs[i],
            &mut cols[i],
        )?);
    }

    // One column buffer for the whole scan (per-sample triplets share
    // `ref_lens[cur_tid]`), the sweep reference, whether it can emit anything
    // at all, and the `-a` bookkeeping fields.
    let mut cur_tid: Option<usize> = None;
    let mut cur_open = false;
    let mut upto = 0u64;
    let mut span_end = 0u64;

    while let Some((mtid, mpos)) = heads.iter().filter_map(|h| h.as_ref()).min().copied() {
        if cur_tid != Some(mtid) {
            // The previous reference still owes trailing zero rows under `-a`.
            if all && cur_open {
                let name = ref_lens[cur_tid.unwrap()].0.as_str();
                zero_fill(
                    &mut out,
                    &mut refs,
                    &mut regions,
                    name,
                    upto + 1,
                    span_end,
                    n,
                )?;
            }
            cur_tid = Some(mtid);
            let (name, len) = &ref_lens[mtid];
            let Some((start, end)) = fill_span(&regions, name, *len) else {
                // No region names this reference: it stays silent, even
                // under `-a`.
                cur_open = false;
                continue;
            };
            cur_open = true;
            upto = start.saturating_sub(1);
            span_end = end;
        }
        let name = ref_lens[mtid].0.as_str();
        let key = (mtid, mpos);

        if all && cur_open {
            // Zero-coverage positions between the last row and this column.
            zero_fill(
                &mut out,
                &mut refs,
                &mut regions,
                name,
                upto + 1,
                (mpos as u64).saturating_sub(1),
                n,
            )?;
        }
        if cur_open && region_allows(&mut regions, name, mpos as u64) {
            let ref_base = match &mut refs {
                Some(cache) => cache
                    .sequence(name)?
                    .and_then(|s| s.get(mpos.saturating_sub(1)))
                    .copied()
                    .unwrap_or(b'N'),
                None => b'N',
            };
            write!(out, "{name}\t{mpos}\t{}", ref_base as char)?;
            // Every input standing on this position contributes its own column;
            // inputs elsewhere are zero triplets, never re-ordered or merged.
            for (col, head) in cols.iter().zip(heads.iter()) {
                if *head == Some(key) {
                    write_sample(&mut out, ref_base, col, args.max_depth as usize)?;
                } else {
                    write!(out, "\t0\t*\t*")?;
                }
            }
            writeln!(out)?;
        }
        upto = mpos as u64;

        // Advance every input standing on this position — whether or not the
        // row was emitted — or the sweep stalls on that input.
        for (((head, col), it), (map, input)) in heads
            .iter_mut()
            .zip(cols.iter_mut())
            .zip(iters.iter_mut())
            .zip(tid_maps.iter().zip(inputs))
        {
            if *head == Some(key) {
                *head = next_col(it, map, input, col)?;
            }
        }
    }

    if all && cur_open {
        let name = ref_lens[cur_tid.unwrap()].0.as_str();
        zero_fill(
            &mut out,
            &mut refs,
            &mut regions,
            name,
            upto + 1,
            span_end,
            n,
        )?;
    }
    out.flush()?;
    Ok(())
}

/// The next column of one input, with its local tid translated into the union
/// tid space. `None` once that input is drained.
fn next_col(
    it: &mut PileupIter,
    map: &[Option<usize>],
    input: &Path,
    col: &mut Vec<(u8, u8)>,
) -> Result<Option<(usize, usize)>> {
    it.next_column_into(col)?
        .map(|(tid, pos)| {
            let utid = map.get(tid).copied().flatten().with_context(|| {
                format!(
                    "{}: reference id {tid} not in the union header",
                    input.display()
                )
            });
            utid.map(|u| (u, pos))
        })
        .transpose()
}

/// One sample's `depth bases quals` triplet with ref-base style markup
/// (`.` for match), capped at the per-sample depth limit.
fn write_sample(
    out: &mut Box<dyn Write>,
    ref_base: u8,
    col: &[(u8, u8)],
    max_depth: usize,
) -> Result<()> {
    let col = &col[..col.len().min(max_depth)];
    write!(out, "\t{}", col.len())?;
    write!(out, "\t")?;
    for &(b, _) in col {
        if b.eq_ignore_ascii_case(&ref_base) {
            out.write_all(b".")?;
        } else {
            out.write_all(&[b.to_ascii_uppercase()])?;
        }
    }
    // Quality column; saturating because the 0xFF "quality unavailable"
    // sentinel overflows a plain `+ 33`.
    write!(out, "\t")?;
    for &(_, q) in col {
        out.write_all(&[q.saturating_add(33).min(126)])?;
    }
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
/// positions. A `0 * *` triplet per input keeps the 3+3N-column mpileup shape,
/// matching samtools' zero-coverage rows. Never emits already-emitted positions.
fn zero_fill(
    out: &mut Box<dyn Write>,
    refs: &mut Option<RefCache>,
    regions: &mut Option<CursorMap>,
    name: &str,
    from: u64,
    to: u64,
    n: usize,
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
        write!(out, "{name}\t{pos}\t{}", ref_base as char)?;
        for _ in 0..n {
            write!(out, "\t0\t*\t*")?;
        }
        writeln!(out)?;
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
        use noodles_sam::header::record::value::map::header::{self as hd, tag};
        // `kira-bam index` refuses headers without `@HD SO:coordinate`, and
        // the indexed-sweep tests build their `.bai` fixtures through it.
        // A builder-built header carries no `@HD` at all, so create it.
        let hd_map = noodles_sam::header::record::value::map::Map::<hd::Header>::builder()
            .insert(tag::SORT_ORDER, hd::sort_order::COORDINATE)
            .build()
            .unwrap();
        let len = NonZeroUsize::try_from(1_000).unwrap();
        Header::builder()
            .set_header(hd_map)
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

    fn write_bam(path: &std::path::Path, header: Header, records: &[RecordBuf]) {
        let mut w = BamWriter::create(Some(path), header, OutputFormat::Bam).unwrap();
        w.write_header().unwrap();
        for r in records {
            w.write_record_buf(r).unwrap();
        }
        w.finish().unwrap();
    }

    /// Write each BAM (with the matching header) — and an `.bai` when `index`
    /// is set — then run mpileup over the list with the given `-r`/`-a`
    /// options, returning every row's tab-separated fields or the run error.
    ///
    /// `index = false` is the classic full-file sweep; `index = true` lets the
    /// eligible-BAM random access kick in, so the two must agree row for row.
    fn rows_multi_opt(
        records: &[&[RecordBuf]],
        headers: &[Header],
        regions: &[&str],
        all: bool,
        min_baseq: u8,
        index: bool,
    ) -> Result<Vec<Vec<String>>> {
        let dir = tempfile::tempdir().unwrap();
        let inputs: Vec<std::path::PathBuf> = records
            .iter()
            .enumerate()
            .map(|(i, rs)| {
                let p = dir.path().join(format!("in{i}.bam"));
                write_bam(p.as_path(), headers[i].clone(), rs);
                if index {
                    crate::index::run(crate::cli::IndexArgs {
                        input: p.clone(),
                        output: None,
                        bai: true,
                        csi: false,
                        min_shift: 14,
                        depth: 5,
                        threads: 1,
                    })
                    .unwrap();
                }
                p
            })
            .collect();
        let out = dir.path().join("out.tsv");
        let res = run(MpileupArgs {
            inputs,
            output: Some(out.clone()),
            reference: None,
            regions: regions.iter().map(|s| s.to_string()).collect(),
            min_mapq: 0,
            min_baseq,
            max_depth: 8000,
            all_positions: all,
        });
        res.map(|_| {
            std::fs::read_to_string(out)
                .unwrap()
                .lines()
                .map(|l| l.split('\t').map(String::from).collect())
                .collect()
        })
    }

    fn rows_multi(
        records: &[&[RecordBuf]],
        headers: &[Header],
        regions: &[&str],
        all: bool,
        min_baseq: u8,
    ) -> Result<Vec<Vec<String>>> {
        rows_multi_opt(records, headers, regions, all, min_baseq, false)
    }

    /// Single input, using the default `chr1`/`chr2` header.
    fn rows(records: &[RecordBuf], regions: &[&str], all: bool) -> Vec<Vec<String>> {
        rows_multi(&[records], &[header()], regions, all, 0).unwrap()
    }

    #[test]
    fn region_restricts_covered_rows() {
        let rows = rows(&[rec("r1", 0, 10)], &["chr1:8-20"], false);
        let pos: Vec<u64> = rows.iter().map(|r| r[1].parse().unwrap()).collect();
        assert_eq!(pos, vec![10, 11, 12, 13, 14]);
        assert!(
            rows.iter().all(|r| r[0] == "chr1" && r[3] == "1"),
            "{rows:?}"
        );
    }

    #[test]
    fn all_positions_fills_region_gaps_and_unions_overlaps() {
        // The union of both regions is discharged once, in order: 1-5 zero,
        // 8-9 zero, then the 5 covered bases, then 15-20 zero.
        let rows = rows(&[rec("r1", 0, 10)], &["chr1:1-5", "chr1:8-20"], true);
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

    /// a covers 10..=14, b covers 12..=16.
    fn two_samples() -> Vec<Vec<String>> {
        rows_multi(
            &[&[rec("a", 0, 10)], &[rec("b", 0, 12)]],
            &[header(), header()],
            &[],
            false,
            0,
        )
        .unwrap()
    }

    #[test]
    fn multiple_inputs_keep_independent_sample_columns() {
        let rows = two_samples();
        // Rows are the union of the two inputs' covered positions; each row
        // carries CHROM POS REF + one depth/bases/quals triplet per input.
        assert_eq!(rows.len(), 7, "{rows:?}");
        assert!(rows.iter().all(|r| r.len() == 9), "{rows:?}");
        let pos: Vec<u64> = rows.iter().map(|r| r[1].parse().unwrap()).collect();
        assert_eq!(pos, vec![10, 11, 12, 13, 14, 15, 16]);
        let depth = |r: &[String], s: usize| r[3 + 3 * s].parse::<u32>().unwrap();
        // Only a covers 10-11; only b covers 15-16; both 12-14.
        assert_eq!(depth(&rows[0], 0), 1);
        assert_eq!(depth(&rows[0], 1), 0);
        assert_eq!(depth(&rows[4], 0), 1);
        assert_eq!(depth(&rows[4], 1), 1);
        assert_eq!(depth(&rows[5], 0), 0);
        assert_eq!(depth(&rows[5], 1), 1);
        // A missing sample is a `0 * *` triplet, not an absent column.
        // Row 0: input b misses (triplet at [6..9]); row 5: input a misses
        // (triplet at [3..6]).
        assert_eq!(&rows[0][6..9], &["0", "*", "*"], "{rows:?}");
        assert_eq!(&rows[5][3..6], &["0", "*", "*"], "{rows:?}");
    }

    #[test]
    fn all_positions_zero_fills_every_sample() {
        // With `-a` the zero rows too must carry one triplet per input.
        let rows = rows_multi(
            &[&[rec("a", 0, 10)], &[rec("b", 0, 12)]],
            &[header(), header()],
            &[],
            true,
            0,
        )
        .unwrap();
        assert_eq!(rows.len(), 1_000, "{rows:?}");
        assert!(rows.iter().all(|r| r.len() == 9), "{rows:?}");
        let depth = |r: &[String], s: usize| r[3 + 3 * s].parse::<u32>().unwrap();
        assert_eq!(depth(&rows[0], 0), 0); // pos 1
        assert_eq!(depth(&rows[9], 0), 1); // pos 10
        assert_eq!(depth(&rows[9], 1), 0);
        assert_eq!(depth(&rows[11], 1), 1); // pos 12
        assert_eq!(depth(&rows[15], 1), 1); // pos 16
        assert_eq!(depth(&rows[15], 0), 0);
    }

    #[test]
    fn sample_columns_follow_cli_input_order() {
        // Inputs reversed: sample 0 is now the 12..=16 read, sample 1 the
        // 10..=14 one, so the depth columns swap exactly.
        let rows = rows_multi(
            &[&[rec("b", 0, 12)], &[rec("a", 0, 10)]],
            &[header(), header()],
            &[],
            false,
            0,
        )
        .unwrap();
        let depth = |r: &[String], s: usize| r[3 + 3 * s].parse::<u32>().unwrap();
        assert_eq!(depth(&rows[0], 0), 0, "{rows:?}"); // pos 10
        assert_eq!(depth(&rows[0], 1), 1);
        assert_eq!(depth(&rows[6], 0), 1, "{rows:?}"); // pos 16
        assert_eq!(depth(&rows[6], 1), 0);
    }

    #[test]
    fn dictionary_union_extends_across_inputs() {
        // The first input knows only chr1; the second adds chr2. chr1 rows
        // come from both samples, chr2 rows only from the second, with a zero
        // triplet for the first sample.
        let h1 = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1_000).unwrap()),
            )
            .build();
        let rows = rows_multi(
            &[&[rec("a", 0, 10)], &[rec("b", 0, 12), rec("c", 1, 12)]],
            &[h1, header()],
            &[],
            false,
            0,
        )
        .unwrap();
        assert!(rows.iter().all(|r| r.len() == 9), "{rows:?}");
        // chr1 portion, then chr2 portion — first dictionary order preserved.
        let chr2_rows: Vec<&Vec<String>> = rows.iter().filter(|r| r[0] == "chr2").collect();
        assert_eq!(chr2_rows.len(), 5, "{rows:?}");
        for r in chr2_rows {
            assert_eq!(&r[3..6], &["0", "*", "*"], "{rows:?}"); // sample a absent
            assert_eq!(r[6], "1"); // sample b covers
        }
    }

    #[test]
    fn incompatible_reordered_dictionary_is_rejected() {
        // A permuted @SQ list cannot be swept in union order; it must error
        // like `merge` does instead of mis-ordering columns.
        let h2 = Header::builder()
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1_000).unwrap()),
            )
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1_000).unwrap()),
            )
            .build();
        let err = rows_multi(
            &[&[rec("a", 0, 10)], &[rec("b", 0, 12)]],
            &[header(), h2],
            &[],
            false,
            0,
        )
        .unwrap_err();
        assert!(
            err.to_string()
                .contains("incompatible with the merged dictionary"),
            "unexpected error: {err:#}"
        );
    }

    /// A 5-base read with 1 bp deletion at its third reference base:
    /// `2M 1D 3M` covering 10..=15.
    fn rec_del(name: &str, start: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Match, 3),
            ]))
            .set_sequence(Sequence::from(b"AAAAA".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 5]))
            .build()
    }

    /// A 5-base perfect match with explicit per-base qualities.
    fn recq(name: &str, start: usize, quals: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, quals.len())]))
            .set_sequence(Sequence::from(vec![b'A'; quals.len()]))
            .set_quality_scores(QualityScores::from(quals.to_vec()))
            .build()
    }

    #[test]
    fn deletion_row_counts_depth_and_prints_the_star() {
        let rows = rows(&[rec_del("r1", 10)], &[], false);
        assert_eq!(
            rows.len(),
            6,
            "the deleted base is an observation: {rows:?}"
        );
        let del = rows.iter().find(|r| r[1] == "12").unwrap();
        assert_eq!(del[0], "chr1");
        assert_eq!(del[3], "1", "deletion contributes to depth: {rows:?}");
        assert_eq!(del[4], "*", "deletion is printed as `*`: {rows:?}");
        assert_eq!(del[5], "!", "deletion carries quality 0: {rows:?}");
    }

    #[test]
    fn base_below_min_baseq_does_not_inflate_depth() {
        let mut quals = vec![30u8; 5];
        quals[0] = 5;
        let rows = rows_multi(&[&[recq("r1", 10, &quals)]], &[header()], &[], false, 13).unwrap();
        assert_eq!(
            rows.len(),
            4,
            "pos 10 is unobserved, not depth 1 with a printed N: {rows:?}",
        );
        assert_eq!(rows[0][1], "11");
        assert!(rows.iter().all(|r| r[3] == "1"), "{rows:?}");
    }

    /// A read jumping a reference skip: `2M 3N 2M` covering 10..=16.
    fn rec_skip(name: &str, start: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Skip, 3),
                Op::new(Kind::Match, 2),
            ]))
            .set_sequence(Sequence::from(b"AAAA".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 4]))
            .build()
    }

    #[test]
    fn skipped_span_is_absent_from_depth_and_output() {
        // 2M covers 10-11, the 3N skip spans 12-14, 2M covers 15-16. The
        // skipped positions produce no rows without `-a` and contribute no
        // deletion observations, unlike a D deletion.
        let rows = rows(&[rec_skip("r1", 10)], &[], false);
        let pos: Vec<u64> = rows.iter().map(|r| r[1].parse().unwrap()).collect();
        assert_eq!(pos, vec![10, 11, 15, 16], "{rows:?}");
        assert!(rows.iter().all(|r| r[3] == "1"), "{rows:?}");
        assert!(
            !rows.iter().any(|r| r[4] == "*"),
            "a skipped position must not print as a deletion: {rows:?}",
        );
    }

    /// Records written as one BAM without an index run through the streaming
    /// sweep; the identical file with an `.bai` runs through index random
    /// access. With `-r` regions the two must produce identical rows — that is
    /// what "indexed instead of streamed, nothing semantic changes" means.
    #[test]
    fn indexed_output_matches_the_streaming_sweep() {
        let records: Vec<RecordBuf> = vec![
            rec("a1", 0, 10),
            rec("a2", 0, 900),
            rec("a3", 1, 500),
            rec("a4", 1, 8),
        ];
        let regions = ["chr1:1-5", "chr1:500-505", "chr2:496-502"];
        for all in [false, true] {
            let streamed = rows_multi(
                &[&records[..3], &records[3..]],
                &[header(), header()],
                &regions,
                all,
                0,
            )
            .unwrap();
            let indexed = rows_multi_opt(
                &[&records[..3], &records[3..]],
                &[header(), header()],
                &regions,
                all,
                0,
                true,
            )
            .unwrap();
            assert_eq!(
                streamed, indexed,
                "indexed rows diverged from the streaming sweep (all={all})"
            );
            assert!(!indexed.is_empty(), "no rows for {regions:?}");
        }
    }

    /// The index must serve a read that starts before the region and crosses
    /// the boundary — only the in-region base is printed, once.
    #[test]
    fn read_spanning_a_region_edge_is_fetched_by_the_index() {
        // `r1` covers 2..=6; region starts at 6, so just position 6 belongs.
        let rows = rows_multi_opt(
            &[&[rec("r1", 0, 2)]],
            &[header()],
            &["chr1:6-20"],
            false,
            0,
            true,
        )
        .unwrap();
        let pos: Vec<u64> = rows.iter().map(|r| r[1].parse().unwrap()).collect();
        assert_eq!(pos, vec![6], "{rows:?}");
        assert!(rows.iter().all(|r| r[3] == "1"), "{rows:?}");
    }

    /// One read covers two disjoint `-r` intervals of the same contig; the
    /// sweep visits its positions in order and the region cursor admits both
    /// intervals' own positions only, so each position is one row with one
    /// observation of this read — never a doubled column.
    #[test]
    fn read_spanning_disjoint_intervals_is_counted_once() {
        // `r1` covers 10..=14: A=[10,11] and B=[14,16]; the read ends exactly
        // at 14 and reaches both, while the covered gap 12..13 has no region.
        let rows = rows_multi_opt(
            &[&[rec("r1", 0, 10)]],
            &[header()],
            &["chr1:10-11", "chr1:14-16"],
            false,
            0,
            true,
        )
        .unwrap();
        let pos: Vec<u64> = rows.iter().map(|r| r[1].parse().unwrap()).collect();
        assert_eq!(pos, vec![10, 11, 14], "{rows:?}");
        assert!(rows.iter().all(|r| r[3] == "1"), "{rows:?}");
    }

    /// A `.bai` that exists but cannot be parsed is an error, never a quiet
    /// fallback to a full scan that would silently change the run.
    #[test]
    fn corrupt_index_reports_an_actionable_error() {
        let dir = tempfile::tempdir().unwrap();
        let bam_path = dir.path().join("in.bam");
        write_bam(&bam_path, header(), &[rec("r1", 0, 10)]);
        std::fs::write(dir.path().join("in.bam.bai"), b"not a bai file").unwrap();

        let err = run(MpileupArgs {
            inputs: vec![bam_path],
            output: Some(dir.path().join("out.tsv")),
            reference: None,
            regions: vec!["chr1:5-15".to_string()],
            min_mapq: 0,
            min_baseq: 0,
            max_depth: 8000,
            all_positions: false,
        })
        .unwrap_err();
        let chain = format!("{err:#}");
        assert!(
            chain.contains("corrupt or stale"),
            "unexpected error: {chain}"
        );
        assert!(
            chain.contains("kira-bam index"),
            "error must suggest the rebuild command: {chain}"
        );
    }

    /// No index at all: `-r` keeps the full-file sweep. Run the same fixture
    /// twice with no `.bai` alongside the input, and once with an index — the
    /// rows must not differ, and the plain BAM must run without an error even
    /// under `-r`.
    #[test]
    fn missing_index_falls_back_to_the_streaming_sweep() {
        let regions = ["chr1:5-15"];
        let records = [rec("r1", 0, 10), rec("r2", 1, 400)];
        // The streaming variant of `rows_multi` writes no index files, so this
        // run exercises the fallback in `run` itself.
        let streamed = rows_multi(&[&records], &[header()], &regions, true, 0).unwrap();
        let indexed = rows_multi_opt(&[&records], &[header()], &regions, true, 0, true).unwrap();
        assert_eq!(streamed, indexed, "identical fixture, divergent rows");
    }

    /// Per-contig cursor regression: records spaced so each sits in its own
    /// binning leaf. `optimize_chunks` merges touching chunk ranges, and in a
    /// coordinate-sorted BAM consecutive records are physically contiguous —
    /// so each contig's index query yields one merged chunk tile, and the
    /// traversal boundary lives between contigs: an `-r` sweep over two
    /// contigs runs one query per contig and seeks between the tiles. Index
    /// geometry is asserted through the same `BinningIndex::query` the
    /// sweep's `IndexedInput` performs, and indexed `mpileup -r` output must
    /// equal the streaming sweep's row for row.
    #[test]
    fn two_contig_index_queries_match_the_streaming_sweep() {
        use noodles_core::{Position, region::Interval};

        // 8 five-base reads per contig, 20 kb apart: each starts in its own
        // 16 kb-wide leaf bin, so every record contributes a chunk.
        const N: usize = 8;
        const STRIDE: usize = 20_000;
        let len = NonZeroUsize::try_from(N * STRIDE).unwrap();
        let contig_len = N * STRIDE;

        use noodles_sam::header::record::value::map::header::{self as hd, tag};
        let hd_map = noodles_sam::header::record::value::map::Map::<hd::Header>::builder()
            .insert(tag::SORT_ORDER, hd::sort_order::COORDINATE)
            .build()
            .unwrap();
        let hdr = Header::builder()
            .set_header(hd_map)
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(len))
            .add_reference_sequence("chr2", Map::<ReferenceSequence>::new(len))
            .build();

        let records = (0..2)
            .flat_map(|tid| (0..N).map(move |i| rec(&format!("r{tid}_{i}"), tid, 1 + i * STRIDE)))
            .collect::<Vec<_>>();

        let dir = tempfile::tempdir().unwrap();
        let bam_path = dir.path().join("in0.bam");
        write_bam(&bam_path, hdr.clone(), &records);

        // True streamed baseline first: with no `.bai` beside the input yet,
        // `run` must take the whole-file sweep. The full row list — not just
        // the positions — is kept for the equality below.
        let regions = vec![
            "chr1:140001-140005".to_string(),
            "chr2:140001-140005".to_string(),
        ];
        let rows_of = |out_path: &std::path::Path| -> Vec<Vec<String>> {
            run(MpileupArgs {
                inputs: vec![bam_path.clone()],
                output: Some(out_path.to_path_buf()),
                reference: None,
                regions: regions.clone(),
                min_mapq: 0,
                min_baseq: 0,
                max_depth: 8000,
                all_positions: false,
            })
            .unwrap();
            std::fs::read_to_string(out_path)
                .unwrap()
                .lines()
                .map(|l| l.split('\t').map(String::from).collect())
                .collect()
        };
        let streamed = rows_of(&dir.path().join("out-stream.tsv"));

        // Assert the index geometry through the same direct query the sweep's
        // `IndexedInput` performs: one merged chunk tile per contig, disjoint
        // from the next, so a two-contig sweep performs one query (and a
        // seek) per contig. Also pin that the spaced records really landed in
        // distinct bins. The index is created only now, after the baseline.
        use noodles_bam as bam;
        use noodles_csi::BinningIndex as _;
        crate::index::run(crate::cli::IndexArgs {
            input: bam_path.clone(),
            output: None,
            bai: true,
            csi: false,
            min_shift: 14,
            depth: 5,
            threads: 1,
        })
        .unwrap();
        let index = bam::bai::fs::read(dir.path().join("in0.bam.bai")).unwrap();
        let whole_contig = Interval::from(
            Position::try_from(1).unwrap()..=Position::try_from(contig_len).unwrap(),
        );
        for tid in 0..2 {
            assert!(
                index.reference_sequences()[tid].bins().len() >= N,
                "the spaced records must occupy distinct bins: {N}"
            );
        }
        let tiled: Vec<_> = (0..2)
            .map(|tid| index.query(tid, whole_contig.clone()).unwrap())
            .collect();
        for (tid, chunks) in tiled.iter().enumerate() {
            assert!(
                !chunks.is_empty(),
                "contig {tid} must hold at least one chunk tile"
            );
        }
        assert!(
            tiled[0]
                .iter()
                .all(|a| tiled[1].iter().all(|b| a.end() <= b.start())),
            "contig chunk tiles must be disjoint: {tiled:?}"
        );
        let total_tiles: usize = tiled.iter().map(|c| c.len()).sum();
        assert!(
            total_tiles >= 2,
            "a two-contig sweep performs one query per contig, got {total_tiles} tiles total"
        );

        // Indexed run now that the `.bai` exists.
        let indexed = rows_of(&dir.path().join("out-indexed.tsv"));
        assert_eq!(
            streamed, indexed,
            "indexed rows diverged from the streaming sweep across \
            the per-contig queries"
        );
        let mut expected = Vec::new();
        for _tid in 0..2 {
            for b in 0..5 {
                expected.push((1 + (N - 1) * STRIDE + b) as u64);
            }
        }
        assert_eq!(
            indexed
                .iter()
                .map(|r| r[1].parse::<u64>().unwrap())
                .collect::<Vec<_>>(),
            expected,
            "{indexed:?}"
        );
        assert!(
            indexed.iter().any(|r| r[0] == "chr1") && indexed.iter().any(|r| r[0] == "chr2"),
            "rows must come from both contigs' queries: {indexed:?}"
        );
    }

    /// BGZF-block regression: the fixture is sized past the 64 KiB cap of one
    /// BGZF block's compressed size (the file-size assert below proves ≥ 2
    /// blocks without any special accessor), so the indexed cursor over the
    /// contig's chunk tile must decode records across block boundaries.
    /// Indexed `mpileup -r` output must still equal the streaming sweep's.
    #[test]
    fn indexed_sweep_reads_across_bgzf_block_boundaries() {
        use noodles_core::Position;

        // 3000 100-base reads, 4 kb apart. Deterministic xorshift bases,
        // qualities and unique names keep deflate from shrinking the fixture.
        const N: usize = 3000;
        const STRIDE: usize = 4096;
        let len = NonZeroUsize::try_from(N * STRIDE).unwrap();

        use noodles_sam::header::record::value::map::header::{self as hd, tag};
        let hd_map = noodles_sam::header::record::value::map::Map::<hd::Header>::builder()
            .insert(tag::SORT_ORDER, hd::sort_order::COORDINATE)
            .build()
            .unwrap();
        let hdr = Header::builder()
            .set_header(hd_map)
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(len))
            .build();

        let records = (0..N)
            .map(|i| {
                let mut state = (i as u64 + 11).wrapping_mul(0x9E37_79B9_7F03_15);
                let mut next = move || {
                    state ^= state >> 12;
                    state ^= state << 25;
                    state ^= state >> 27;
                    state.wrapping_mul(0x2545_F491_4F6C_DD1D)
                };
                let seq: Vec<u8> = (0..100).map(|_| b"ACGT"[(next() & 3) as usize]).collect();
                let quals: Vec<u8> = (0..100).map(|_| ((next() & 0x3f) + 30) as u8).collect();
                RecordBuf::builder()
                    .set_name(format!("q{i}"))
                    .set_flags(Flags::default())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(1 + i * STRIDE).unwrap())
                    .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 100)]))
                    .set_sequence(Sequence::from(seq))
                    .set_quality_scores(QualityScores::from(quals))
                    .build()
            })
            .collect::<Vec<_>>();

        let dir = tempfile::tempdir().unwrap();
        let bam_path = dir.path().join("in0.bam");
        write_bam(&bam_path, hdr.clone(), &records);

        // One BGZF block's `BSIZE` — its compressed footprint, header/body
        // included — is capped at 65536 bytes, so a larger file guarantees ≥ 2
        // blocks. The high-entropy qualities (3000 * 100 random bytes) make
        // the compressed size comfortably exceed the cap.
        let file_size = std::fs::metadata(&bam_path).unwrap().len();
        assert!(
            file_size > 64 * 1024,
            "fixture must hold >1 BGZF block, got {file_size} bytes"
        );

        // True streamed baseline first: with no `.bai` beside the input yet,
        // `run` must take the whole-file sweep. Full rows are kept for the
        // equality below. The index is created only after this call.
        let regions = vec!["chr1:4096001-4100196".to_string()];
        let rows_of = |out_path: &std::path::Path| -> Vec<Vec<String>> {
            run(MpileupArgs {
                inputs: vec![bam_path.clone()],
                output: Some(out_path.to_path_buf()),
                reference: None,
                regions: regions.clone(),
                min_mapq: 0,
                min_baseq: 0,
                max_depth: 8000,
                all_positions: false,
            })
            .unwrap();
            std::fs::read_to_string(out_path)
                .unwrap()
                .lines()
                .map(|l| l.split('\t').map(String::from).collect())
                .collect()
        };
        let streamed = rows_of(&dir.path().join("out-stream.tsv"));

        // Index build + same direct query sanity as the sibling tests.
        use noodles_bam as bam;
        use noodles_csi::BinningIndex as _;
        crate::index::run(crate::cli::IndexArgs {
            input: bam_path.clone(),
            output: None,
            bai: true,
            csi: false,
            min_shift: 14,
            depth: 5,
            threads: 1,
        })
        .unwrap();
        let index = bam::bai::fs::read(dir.path().join("in0.bam.bai")).unwrap();
        let whole_contig = noodles_core::region::Interval::from(
            Position::try_from(1).unwrap()..=Position::try_from(N * STRIDE).unwrap(),
        );
        let tiles = index.query(0, whole_contig).unwrap();
        assert!(
            !tiles.is_empty(),
            "the contig's chunk tile list must not be empty"
        );

        // Indexed run now that the `.bai` exists.
        let indexed = rows_of(&dir.path().join("out-indexed.tsv"));
        assert_eq!(
            streamed, indexed,
            "indexed rows diverged from the streaming sweep across BGZF blocks"
        );
        assert_eq!(indexed.len(), 200, "{indexed:?}");
        assert_eq!(indexed[0][1], "4096001");
        assert_eq!(indexed[100][1], "4100097");
        assert_eq!(indexed[199][1], "4100196");
    }
}
