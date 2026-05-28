use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_csi::{self as csi, BinningIndex};

use crate::cli::IdxstatsArgs;

pub fn run(args: IdxstatsArgs) -> Result<()> {
    let bam_path = &args.input;
    let mut reader = bam::io::Reader::new(File::open(bam_path).context("open BAM")?);
    let header = reader.read_header().context("read header")?;

    // Try .bai first, fall back to .csi.
    let index: Box<dyn BinningIndex> = {
        let bai_path = build_index_path(bam_path, "bai");
        let csi_path = build_index_path(bam_path, "csi");
        if bai_path.exists() {
            Box::new(bam::bai::fs::read(bai_path).context("read .bai")?)
        } else if csi_path.exists() {
            Box::new(csi::fs::read(csi_path).context("read .csi")?)
        } else {
            anyhow::bail!("no .bai or .csi index found for {}", bam_path.display());
        }
    };

    let mut out: Box<dyn Write> = Box::new(BufWriter::new(std::io::stdout()));

    // Each reference: name \t length \t mapped \t unmapped
    let refs: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .map(|(name, sq)| (name.to_string(), usize::from(sq.length())))
        .collect();

    let mut ref_iter = index.reference_sequences();
    for (name, len) in &refs {
        let (mapped, unmapped) = ref_iter
            .next()
            .and_then(|r| r.metadata())
            .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
            .unwrap_or((0, 0));
        writeln!(out, "{name}\t{len}\t{mapped}\t{unmapped}")?;
    }

    // Unplaced (samtools row: `*\t0\t0\tN`)
    let unplaced = index.unplaced_unmapped_record_count().unwrap_or(0);
    writeln!(out, "*\t0\t0\t{unplaced}")?;
    out.flush()?;
    Ok(())
}

fn build_index_path(bam: &std::path::Path, ext: &str) -> std::path::PathBuf {
    let mut s = bam.as_os_str().to_os_string();
    s.push(".");
    s.push(ext);
    std::path::PathBuf::from(s)
}
