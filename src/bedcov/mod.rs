use std::collections::HashSet;
use std::fs::File;
use std::io::Write;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::{Position, Region};
use noodles_sam::alignment::Record as _;

use crate::bed::{self, Interval};
use crate::cli::BedcovArgs;
use crate::io::BufferedOutput;

pub fn run(args: BedcovArgs) -> Result<()> {
    let intervals = bed::read(&args.bed)?;
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
            if u16::from(r.flags()) & 0x704 != 0 {
                continue;
            }
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
