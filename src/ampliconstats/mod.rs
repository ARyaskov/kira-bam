use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Region;
use noodles_sam::alignment::RecordBuf;

use crate::cli::AmpliconstatsArgs;

#[derive(Debug, Clone)]
struct Amp {
    name: String,
    chrom: String,
    start: u64,
    end: u64,
}

pub fn run(args: AmpliconstatsArgs) -> Result<()> {
    let amps = read_amplicons(&args.bed)?;

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    writeln!(out, "# amplicon\tchrom\tstart\tend\tlength\treads...")?;
    for amp in &amps {
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            amp.name,
            amp.chrom,
            amp.start,
            amp.end,
            amp.end.saturating_sub(amp.start) + 1
        )?;
        for input in &args.inputs {
            let n = count_in_region(input, amp, args.min_mapq).unwrap_or(0);
            write!(out, "\t{n}")?;
        }
        writeln!(out)?;
    }
    Ok(())
}

fn read_amplicons(path: &std::path::Path) -> Result<Vec<Amp>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out = Vec::new();
    let mut idx = 0;
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        idx += 1;
        let name = if parts.len() >= 4 {
            parts[3].to_string()
        } else {
            format!("amp{idx:04}")
        };
        out.push(Amp {
            name,
            chrom: parts[0].to_string(),
            start: parts[1].parse::<u64>().unwrap_or(0) + 1,
            end: parts[2].parse().unwrap_or(0),
        });
    }
    Ok(out)
}

fn count_in_region(input: &std::path::Path, amp: &Amp, min_mapq: u8) -> Result<u64> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(input)
        .context("open indexed BAM")?;
    let header = reader.read_header()?;
    let region: Region = format!("{}:{}-{}", amp.chrom, amp.start, amp.end)
        .parse()
        .map_err(|e| anyhow::anyhow!("region parse: {e:?}"))?;
    let q = reader.query(&header, &region)?;
    let mut n = 0u64;
    for r in q.records() {
        let r = r?;
        let buf = RecordBuf::try_from_alignment_record(&header, &r)?;
        if let Some(mq) = buf.mapping_quality()
            && u8::from(mq) < min_mapq
        {
            continue;
        }
        n += 1;
    }
    Ok(n)
}
