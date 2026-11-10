use std::fs::File;
use std::io::{BufWriter, Write};

use std::collections::HashMap;
use std::io::BufRead;

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use rustc_hash::FxHasher;

use crate::cli::StatsArgs;
use crate::io::{BamReader, install_thread_pool};

#[derive(Debug, Clone)]
struct BedIv {
    start: u64,
    end: u64,
}

fn load_bed(
    path: &std::path::Path,
) -> Result<HashMap<String, Vec<BedIv>, std::hash::BuildHasherDefault<FxHasher>>> {
    let f = std::fs::File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out: HashMap<String, Vec<BedIv>, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    for line in std::io::BufReader::new(f).lines() {
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
        out.entry(parts[0].to_string()).or_default().push(BedIv {
            start: start + 1,
            end,
        });
    }
    Ok(out)
}

fn in_any(ivs: &[BedIv], start: u64, end: u64) -> bool {
    ivs.iter().any(|iv| iv.start <= end && iv.end >= start)
}

const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAP: u16 = 0x4;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QCFAIL: u16 = 0x200;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Initial histogram capacities, not caps. They used to be hard limits, so any
/// insert or read longer than this vanished from the report with no indication
/// — silently truncating exactly the tail a QC plot is there to show. Long-read
/// data made both routinely wrong.
const INITIAL_INSERT_SIZE_BINS: usize = 8000;
const INITIAL_READ_LENGTH_BINS: usize = 1024;

/// Record a hit in `hist` at `idx`, growing it as needed.
fn bump(hist: &mut Vec<u64>, idx: usize) {
    if idx >= hist.len() {
        hist.resize(idx + 1, 0);
    }
    hist[idx] += 1;
}

#[derive(Default)]
struct Stats {
    raw_total: u64,
    filtered_out: u64,
    mapped: u64,
    duplicates: u64,
    paired: u64,
    properly_paired: u64,
    qc_failed: u64,
    bases_total: u64,
    /// "bases mapped": total read length summed over mapped reads.
    bases_mapped_seq: u64,
    /// "bases mapped (cigar)": only reference-and-read-consuming CIGAR ops.
    bases_mapped: u64,
    bases_n: u64,
    reads_unmapped: u64,
    reads_qc_failed: u64,
    reads_mq0: u64,
    insert_size_hist: Vec<u64>,
    read_length_hist: Vec<u64>,
    gc_hist: Vec<u64>,
    mapq_hist: Vec<u64>,
    sum_mapq: u64,
}

pub fn run(args: StatsArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let bed = match &args.bed {
        Some(p) => Some(load_bed(p)?),
        None => None,
    };

    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let header = reader.header().clone();
    let mut st = Stats {
        insert_size_hist: vec![0; INITIAL_INSERT_SIZE_BINS + 1],
        read_length_hist: vec![0; INITIAL_READ_LENGTH_BINS + 1],
        gc_hist: vec![0; 101],
        mapq_hist: vec![0; 256],
        ..Default::default()
    };
    // Contig names resolved once, not rebuilt as a String per record.
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| n.to_string())
        .collect();

    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        let flags = u16::from(rec.flags());
        st.raw_total += 1;
        if (flags & FLAG_SECONDARY) != 0 || (flags & FLAG_SUPPLEMENTARY) != 0 {
            st.filtered_out += 1;
            continue;
        }
        // BED filter (samtools issue #2172).
        if let Some(b) = &bed {
            if let Some(tid) = rec.reference_sequence_id() {
                let ref_name = ref_names.get(tid).map(String::as_str).unwrap_or("");
                let s = rec
                    .alignment_start()
                    .map(|p| usize::from(p) as u64)
                    .unwrap_or(0);
                let e = rec
                    .alignment_end()
                    .map(|p| usize::from(p) as u64)
                    .unwrap_or(s);
                let pass = b
                    .get(ref_name)
                    .map(|ivs| in_any(ivs, s, e))
                    .unwrap_or(false);
                if !pass {
                    st.filtered_out += 1;
                    continue;
                }
            } else {
                st.filtered_out += 1;
                continue;
            }
        }
        // MAPQ 255 means "unavailable", which noodles reports as `None`. Those
        // records used to skip the `-q` filter and every MAPQ counter entirely;
        // treating the sentinel as its literal 255 keeps one record on one code
        // path.
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(255);
        if mapq < args.min_mapq {
            st.filtered_out += 1;
            continue;
        }
        st.mapq_hist[mapq as usize] += 1;
        if mapq == 0 {
            st.reads_mq0 += 1;
        }
        if (flags & FLAG_QCFAIL) != 0 {
            st.qc_failed += 1;
            st.reads_qc_failed += 1;
        }
        if (flags & FLAG_DUP) != 0 {
            st.duplicates += 1;
        }
        if (flags & FLAG_UNMAP) != 0 {
            st.reads_unmapped += 1;
        } else {
            st.mapped += 1;
            // Summed over the same records "average MAPQ" divides by, so the
            // two agree; unmapped reads used to contribute to the numerator
            // while the denominator counted only mapped ones.
            st.sum_mapq += mapq as u64;
            // "bases mapped" = full read length of mapped reads.
            st.bases_mapped_seq += rec.sequence().as_ref().len() as u64;
            // "bases mapped (cigar)" = only ref-and-read-consuming ops.
            let mapped_bases: u64 = rec
                .cigar()
                .as_ref()
                .iter()
                .filter(|op| op.kind().consumes_reference() && op.kind().consumes_read())
                .map(|op| op.len() as u64)
                .sum();
            st.bases_mapped += mapped_bases;
        }
        if (flags & FLAG_PAIRED) != 0 {
            st.paired += 1;
            if (flags & FLAG_PROPER_PAIR) != 0 {
                st.properly_paired += 1;
            }
            // One count per *pair*, not per read: samtools tallies the insert
            // size only on the record with a positive TLEN, so counting both
            // mates doubled every bin of the IS histogram.
            let tlen = rec.template_length();
            if tlen > 0 {
                bump(&mut st.insert_size_hist, tlen as usize);
            }
        }
        let seq = rec.sequence().as_ref();
        let read_len = seq.len();
        st.bases_total += read_len as u64;
        bump(&mut st.read_length_hist, read_len);
        let n_count = seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
        st.bases_n += n_count as u64;
        let gc = seq
            .iter()
            .filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
            .count();
        if let Some(gc_pct) = (100 * gc).checked_div(read_len) {
            st.gc_hist[gc_pct.min(100)] += 1;
        }
    }

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    if let Some(s) = &args.sample_name {
        writeln!(out, "# sample {s}")?;
    }
    writeln!(out, "# This file was produced by kira-bam stats.")?;
    // Field names and order follow samtools stats, which is what MultiQC and
    // every other `SN`-line parser keys on. `total bases` in particular is
    // spelled `total length` there and was simply not found before.
    writeln!(out, "SN\traw total sequences:\t{}", st.raw_total)?;
    writeln!(out, "SN\tfiltered sequences:\t{}", st.filtered_out)?;
    writeln!(out, "SN\tsequences:\t{}", st.raw_total - st.filtered_out)?;
    writeln!(out, "SN\treads mapped:\t{}", st.mapped)?;
    writeln!(out, "SN\treads unmapped:\t{}", st.reads_unmapped)?;
    writeln!(out, "SN\treads properly paired:\t{}", st.properly_paired)?;
    writeln!(out, "SN\treads paired:\t{}", st.paired)?;
    writeln!(out, "SN\treads duplicated:\t{}", st.duplicates)?;
    writeln!(out, "SN\treads MQ0:\t{}", st.reads_mq0)?;
    writeln!(out, "SN\treads QC failed:\t{}", st.reads_qc_failed)?;
    writeln!(out, "SN\ttotal length:\t{}", st.bases_total)?;
    writeln!(out, "SN\tbases mapped:\t{}", st.bases_mapped_seq)?;
    writeln!(out, "SN\tbases mapped (cigar):\t{}", st.bases_mapped)?;
    writeln!(out, "SN\tbases N:\t{}", st.bases_n)?;
    writeln!(
        out,
        "SN\taverage MAPQ:\t{:.2}",
        if st.mapped > 0 {
            st.sum_mapq as f64 / st.mapped as f64
        } else {
            0.0
        }
    )?;

    for (len, &n) in st.read_length_hist.iter().enumerate() {
        if n > 0 {
            writeln!(out, "RL\t{len}\t{n}")?;
        }
    }
    for (isize, &n) in st.insert_size_hist.iter().enumerate() {
        if n > 0 {
            writeln!(out, "IS\t{isize}\t{n}")?;
        }
    }
    for (gc, &n) in st.gc_hist.iter().enumerate() {
        if n > 0 {
            writeln!(out, "GCF\t{gc}\t{n}")?;
        }
    }
    for (mq, &n) in st.mapq_hist.iter().enumerate() {
        if n > 0 {
            writeln!(out, "MAPQ\t{mq}\t{n}")?;
        }
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::BamWriter;
    use crate::types::OutputFormat;
    use noodles_core::Position;
    use noodles_sam::Header;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles_sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    fn header() -> Header {
        Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1_000).unwrap()),
            )
            .build()
    }

    /// A 5 bp perfect match with explicit flags.
    fn rec(name: &str, flags: Flags) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5)]))
            .set_sequence(Sequence::from(b"ACGTA".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 5]))
            .build()
    }

    /// `SN <label> value` counts from a full run over a BAM of `records`.
    fn sn_counts(records: &[RecordBuf]) -> Vec<(String, String)> {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("in.bam");
        let out = dir.path().join("out.txt");
        {
            let mut w = BamWriter::create(Some(&bam), header(), OutputFormat::Bam).unwrap();
            w.write_header().unwrap();
            for r in records {
                w.write_record_buf(r).unwrap();
            }
            w.finish().unwrap();
        }
        run(StatsArgs {
            input: bam,
            output: Some(out.clone()),
            reference: None,
            bed: None,
            sample_name: None,
            min_mapq: 0,
            threads: 1,
        })
        .unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        text.lines()
            .filter_map(|l| {
                let mut f = l.split('\t');
                match f.next() {
                    Some("SN") => Some((
                        f.next().unwrap().trim_end_matches(':').to_string(),
                        f.next().unwrap().to_string(),
                    )),
                    _ => None,
                }
            })
            .collect()
    }

    #[test]
    fn secondary_and_supplementary_records_are_counted_as_filtered() {
        let counts = sn_counts(&[
            rec("r1", Flags::default()),
            rec("r2", Flags::SECONDARY),
            rec("r3", Flags::SUPPLEMENTARY),
        ]);
        let counter = |label: &str| {
            counts
                .iter()
                .find(|(l, _)| l == label)
                .map(|(_, v)| {
                    v.parse::<u64>()
                        .unwrap_or_else(|_| panic!("non-integer {label}: {v}"))
                })
                .unwrap_or_else(|| panic!("missing {label} line: {counts:?}"))
        };
        // All 3 records are raw; the 2 flagged ones are filtered and the
        // summary `sequences` count only keeps the primary record.
        assert_eq!(counter("raw total sequences"), 3, "{counts:?}");
        assert_eq!(counter("filtered sequences"), 2, "{counts:?}");
        assert_eq!(counter("sequences"), 1, "{counts:?}");
    }

    /// The documented reference contract, run without `-r`: the GC histogram
    /// is unconditional and samtools' error-rate sections (FFQ/LRQ) are never
    /// written — the help and README must keep saying exactly that.
    #[test]
    fn gc_is_emitted_without_a_reference_and_no_error_rate_section_exists() {
        let counts = sn_counts(&[rec("r1", Flags::default())]);
        assert!(
            counts.iter().any(|(l, _)| l.contains("average MAPQ")),
            "the run produced no SN summary: {counts:?}"
        );

        let raw = gc_report_text(&[rec("r1", Flags::default())]);
        assert!(raw.lines().any(|l| l.starts_with("GCF\t")));
        assert!(
            !raw.lines()
                .any(|l| l.starts_with("FFQ\t") || l.starts_with("LRQ\t")),
            "an error-rate section was emitted: {raw:?}"
        );
    }

    /// Full report text (SN/RL/IS/GCF/MAPQ lines) from one run, like
    /// [`sn_counts`] but without the SN filter.
    fn gc_report_text(records: &[RecordBuf]) -> String {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("in.bam");
        let out = dir.path().join("out.txt");
        {
            let mut w = BamWriter::create(Some(&bam), header(), OutputFormat::Bam).unwrap();
            w.write_header().unwrap();
            for r in records {
                w.write_record_buf(r).unwrap();
            }
            w.finish().unwrap();
        }
        run(StatsArgs {
            input: bam,
            output: Some(out.clone()),
            reference: None,
            bed: None,
            sample_name: None,
            min_mapq: 0,
            threads: 1,
        })
        .unwrap();
        std::fs::read_to_string(out).unwrap()
    }

    /// Pin the help text's truthfulness about `-r` — it must match the
    /// behavior above (CRAM-only reference, unconditional GC, no error rate).
    #[test]
    fn stats_help_states_the_reference_contract() {
        use clap::CommandFactory;
        use std::fmt::Write as _;
        let mut cmd = crate::cli::StatsArgs::command();
        let mut text = String::new();
        // `render_help`'s plain write is enough; no wall-of-text snapshot.
        write!(text, "{}", cmd.render_long_help()).unwrap();
        for needle in ["CRAM", "unconditional", "error-rate"] {
            assert!(
                text.contains(needle),
                "stats help must mention {needle:?}; got:\n{text}"
            );
        }
    }
}
