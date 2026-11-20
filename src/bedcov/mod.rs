use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use anyhow::{Context, Result, ensure};
use noodles_bam as bam;
use noodles_core::{Position, Region};
use noodles_sam::alignment::Record as _;

use crate::cli::BedcovArgs;
use crate::io::BufferedOutput;

#[derive(Debug, Clone)]
struct Interval {
    chrom: String,
    start: u64,
    end: u64,
}

pub fn run(args: BedcovArgs) -> Result<()> {
    let intervals = read_bed(&args.bed)?;
    let mut out = BufferedOutput::create(args.output.as_deref())?;

    // Open each input once. Reopening per interval meant a file open, a header
    // parse and a full index load for every cell of the output — quadratic in
    // the work that actually matters on a panel with thousands of targets.
    let mut sources: Vec<Source> = args
        .inputs
        .iter()
        .map(|p| Source::open(p))
        .collect::<Result<_>>()?;

    // Contigs the BED names but a BAM does not have: worth one warning each,
    // not a silent column of zeros repeated per interval.
    let mut warned: HashSet<(usize, String)> = HashSet::new();

    for iv in &intervals {
        write!(out, "{}\t{}\t{}", iv.chrom, iv.start, iv.end)?;
        for (i, src) in sources.iter_mut().enumerate() {
            if !src.has_reference(&iv.chrom) {
                if warned.insert((i, iv.chrom.clone())) {
                    eprintln!(
                        "[kira-bam bedcov] {}: no @SQ named {} — reporting 0",
                        args.inputs[i].display(),
                        iv.chrom
                    );
                }
                write!(out, "\t0")?;
                continue;
            }
            // Real query failures (a corrupt index, an I/O error) now abort
            // instead of being flattened into a zero that reads as a valid
            // measurement.
            let cov = src.count_coverage(iv, args.min_mapq).with_context(|| {
                format!(
                    "count coverage of {}:{}-{} in {}",
                    iv.chrom,
                    iv.start,
                    iv.end,
                    args.inputs[i].display()
                )
            })?;
            write!(out, "\t{cov}")?;
        }
        writeln!(out)?;
    }
    out.finish()
}

/// One input BAM, held open with its header and index for the whole run.
struct Source {
    reader: bam::io::IndexedReader<noodles_bgzf::io::Reader<File>>,
    header: noodles_sam::Header,
}

impl Source {
    fn open(path: &std::path::Path) -> Result<Self> {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("open indexed BAM {}", path.display()))?;
        let header = reader.read_header().context("read header")?;
        Ok(Self { reader, header })
    }

    fn has_reference(&self, name: &str) -> bool {
        self.header
            .reference_sequences()
            .contains_key(name.as_bytes())
    }

    fn count_coverage(&mut self, iv: &Interval, min_mapq: u8) -> Result<u64> {
        if iv.start == iv.end {
            return Ok(0);
        }
        let start = Position::try_from(usize::try_from(iv.start + 1)?)?;
        let end = Position::try_from(usize::try_from(iv.end)?)?;
        let region = Region::new(iv.chrom.as_bytes(), start..=end);
        let query = self.reader.query(&self.header, &region)?;
        let mut sum: u64 = 0;
        for r in query.records() {
            let r = r?;
            // The raw record's fixed fields are all this needs; materialising a
            // whole `RecordBuf` decoded the name, CIGAR, sequence and tags only
            // to drop them.
            if let Some(mapq) = r.mapping_quality()
                && u8::from(mapq) < min_mapq
            {
                continue;
            }
            let (Some(s), Some(e)) = (
                r.alignment_start().transpose()?,
                r.alignment_end().transpose()?,
            ) else {
                continue;
            };
            let s = usize::from(s) as u64;
            let e = usize::from(e) as u64;
            let overlap_start = s.max(iv.start + 1);
            let overlap_end = e.min(iv.end);
            if overlap_end >= overlap_start {
                sum += overlap_end - overlap_start + 1;
            }
        }
        Ok(sum)
    }
}

fn read_bed(path: &std::path::Path) -> Result<Vec<Interval>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    for (index, line) in r.lines().enumerate() {
        let location = || format!("{}: line {}", path.display(), index + 1);
        let line = line.with_context(location)?;
        if line.trim().is_empty()
            || line.starts_with('#')
            || line.starts_with("track ")
            || line.starts_with("browser ")
        {
            continue;
        }
        out.push(parse_interval(&line).with_context(location)?);
    }
    Ok(out)
}

fn parse_interval(line: &str) -> Result<Interval> {
    let mut fields = line.split('\t');
    let chrom = fields
        .next()
        .filter(|field| !field.is_empty())
        .context("missing chromosome")?;
    let start = fields
        .next()
        .context("missing BED start")?
        .parse()
        .context("invalid BED start")?;
    let end = fields
        .next()
        .context("missing BED end")?
        .parse()
        .context("invalid BED end")?;
    ensure!(start <= end, "BED start exceeds end");
    ensure!(
        usize::try_from(end).is_ok(),
        "BED end exceeds addressable coordinates"
    );
    Ok(Interval {
        chrom: chrom.to_owned(),
        start,
        end,
    })
}
