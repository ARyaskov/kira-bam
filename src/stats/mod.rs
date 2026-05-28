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

const MAX_INSERT_SIZE: usize = 8000;
const MAX_READ_LENGTH: usize = 1024;

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

    let mut reader = BamReader::open(&args.input).context("open input")?;
    let header = reader.header().clone();
    let mut st = Stats {
        insert_size_hist: vec![0; MAX_INSERT_SIZE + 1],
        read_length_hist: vec![0; MAX_READ_LENGTH + 1],
        gc_hist: vec![0; 101],
        mapq_hist: vec![0; 256],
        ..Default::default()
    };

    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        let flags = u16::from(rec.flags());
        st.raw_total += 1;
        if (flags & FLAG_SECONDARY) != 0 || (flags & FLAG_SUPPLEMENTARY) != 0 {
            continue;
        }
        // BED filter (samtools issue #2172).
        if let Some(b) = &bed {
            if let Some(tid) = rec.reference_sequence_id() {
                let ref_name: String = header
                    .reference_sequences()
                    .get_index(tid)
                    .map(|(n, _)| n.to_string())
                    .unwrap_or_default();
                let s = rec
                    .alignment_start()
                    .map(|p| usize::from(p) as u64)
                    .unwrap_or(0);
                let e = rec
                    .alignment_end()
                    .map(|p| usize::from(p) as u64)
                    .unwrap_or(s);
                let pass = b
                    .get(&ref_name)
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
        if let Some(mapq) = rec.mapping_quality() {
            let mapq = u8::from(mapq);
            if mapq < args.min_mapq {
                st.filtered_out += 1;
                continue;
            }
            st.sum_mapq += mapq as u64;
            st.mapq_hist[mapq as usize] += 1;
            if mapq == 0 {
                st.reads_mq0 += 1;
            }
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
            // Count mapped bases via CIGAR.
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
            let isize = rec.template_length().unsigned_abs() as usize;
            if isize <= MAX_INSERT_SIZE {
                st.insert_size_hist[isize] += 1;
            }
        }
        let seq = rec.sequence().as_ref();
        let read_len = seq.len();
        st.bases_total += read_len as u64;
        if read_len <= MAX_READ_LENGTH {
            st.read_length_hist[read_len] += 1;
        }
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
    writeln!(out, "SN\traw total sequences:\t{}", st.raw_total)?;
    writeln!(out, "SN\tfiltered sequences:\t{}", st.filtered_out)?;
    writeln!(out, "SN\treads mapped:\t{}", st.mapped)?;
    writeln!(out, "SN\treads unmapped:\t{}", st.reads_unmapped)?;
    writeln!(out, "SN\treads duplicated:\t{}", st.duplicates)?;
    writeln!(out, "SN\treads QC failed:\t{}", st.reads_qc_failed)?;
    writeln!(out, "SN\treads paired:\t{}", st.paired)?;
    writeln!(out, "SN\treads properly paired:\t{}", st.properly_paired)?;
    writeln!(out, "SN\treads MQ0:\t{}", st.reads_mq0)?;
    writeln!(out, "SN\ttotal bases:\t{}", st.bases_total)?;
    writeln!(out, "SN\tbases mapped:\t{}", st.bases_mapped)?;
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
