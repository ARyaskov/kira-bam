use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Region;
use noodles_sam::alignment::RecordBuf;

use crate::cli::BedcovArgs;

#[derive(Debug, Clone)]
struct Interval {
    chrom: String,
    start: u64,
    end: u64,
}

pub fn run(args: BedcovArgs) -> Result<()> {
    let intervals = read_bed(&args.bed)?;
    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    for iv in &intervals {
        write!(out, "{}\t{}\t{}", iv.chrom, iv.start, iv.end)?;
        for bam_path in &args.inputs {
            let cov = count_coverage(bam_path, iv, args.min_mapq).unwrap_or(0);
            write!(out, "\t{cov}")?;
        }
        writeln!(out)?;
    }
    out.flush()?;
    Ok(())
}

fn read_bed(path: &std::path::Path) -> Result<Vec<Interval>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    for line in r.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let start: u64 = parts[1].parse().unwrap_or(0);
        let end: u64 = parts[2].parse().unwrap_or(0);
        out.push(Interval {
            chrom: parts[0].to_string(),
            start,
            end,
        });
    }
    Ok(out)
}

fn count_coverage(bam_path: &std::path::Path, iv: &Interval, min_mapq: u8) -> Result<u64> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .context("open indexed BAM")?;
    let header = reader.read_header()?;
    let region: Region = format!("{}:{}-{}", iv.chrom, iv.start + 1, iv.end)
        .parse()
        .map_err(|e| anyhow::anyhow!("region parse: {e:?}"))?;
    let query = reader.query(&header, &region)?;
    let mut sum: u64 = 0;
    for r in query.records() {
        let r = r?;
        let buf = RecordBuf::try_from_alignment_record(&header, &r)?;
        if let Some(mapq) = buf.mapping_quality()
            && u8::from(mapq) < min_mapq
        {
            continue;
        }
        if let (Some(s), Some(e)) = (buf.alignment_start(), buf.alignment_end()) {
            let s = usize::from(s) as u64;
            let e = usize::from(e) as u64;
            let overlap_start = s.max(iv.start + 1);
            let overlap_end = e.min(iv.end);
            if overlap_end >= overlap_start {
                sum += overlap_end - overlap_start + 1;
            }
        }
    }
    Ok(sum)
}
