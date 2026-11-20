use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Region;
use noodles_sam::Header;

use crate::cli::AmpliconstatsArgs;
use crate::io::BufferedOutput;

#[derive(Debug, Clone)]
struct Amp {
    name: String,
    chrom: String,
    start: u64,
    end: u64,
}

pub fn run(args: AmpliconstatsArgs) -> Result<()> {
    let amps = read_amplicons(&args.bed)?;
    let mut sources = args
        .inputs
        .iter()
        .map(|path| {
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .with_context(|| format!("open indexed BAM {}", path.display()))?;
            let header = reader
                .read_header()
                .with_context(|| format!("read header of {}", path.display()))?;
            Ok((reader, header))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut out = BufferedOutput::create(args.output.as_deref())?;

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
        for (input, (reader, header)) in args.inputs.iter().zip(&mut sources) {
            let n = count_in_region(reader, header, amp, args.min_mapq).with_context(|| {
                format!(
                    "count reads in {}:{}-{} in {}",
                    amp.chrom,
                    amp.start,
                    amp.end,
                    input.display()
                )
            })?;
            write!(out, "\t{n}")?;
        }
        writeln!(out)?;
    }
    out.finish()
}

fn read_amplicons(path: &std::path::Path) -> Result<Vec<Amp>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out = Vec::new();
    let mut idx = 0;
    for (line_number, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track ")
            || line.starts_with("browser ")
        {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            anyhow::bail!(
                "{}:{}: expected three BED columns",
                path.display(),
                line_number + 1
            );
        }
        idx += 1;
        let name = if parts.len() >= 4 {
            parts[3].to_string()
        } else {
            format!("amp{idx:04}")
        };
        let bounds = || {
            format!(
                "{}:{}: invalid BED interval",
                path.display(),
                line_number + 1
            )
        };
        let start = parts[1].parse::<u64>().with_context(bounds)?;
        let end = parts[2].parse::<u64>().with_context(bounds)?;
        anyhow::ensure!(start < end, "{}", bounds());
        out.push(Amp {
            name,
            chrom: parts[0].to_string(),
            start: start + 1,
            end,
        });
    }
    Ok(out)
}

fn count_in_region(
    reader: &mut bam::io::IndexedReader<noodles_bgzf::io::Reader<File>>,
    header: &Header,
    amp: &Amp,
    min_mapq: u8,
) -> Result<u64> {
    let region: Region = format!("{}:{}-{}", amp.chrom, amp.start, amp.end)
        .parse()
        .map_err(|e| anyhow::anyhow!("region parse: {e:?}"))?;
    let q = reader.query(header, &region)?;
    let mut n = 0u64;
    for r in q.records() {
        let r = r?;
        if let Some(mq) = r.mapping_quality()
            && u8::from(mq) < min_mapq
        {
            continue;
        }
        n += 1;
    }
    Ok(n)
}
