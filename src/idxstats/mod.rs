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

    anyhow::ensure!(
        index.reference_sequences().count() == refs.len(),
        "index reference count does not match the BAM header"
    );
    let indexed_counts = index
        .reference_sequences()
        .map(|reference| {
            reference
                .metadata()
                .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
        })
        .collect::<Option<Vec<_>>>();
    let (counts, unplaced) = match (indexed_counts, index.unplaced_unmapped_record_count()) {
        (Some(counts), Some(unplaced)) => (counts, unplaced),
        _ => {
            let mut counts = vec![(0, 0); refs.len()];
            let mut unplaced = 0;
            let mut record = bam::Record::default();
            while reader.read_record(&mut record)? != 0 {
                match record.reference_sequence_id().transpose()? {
                    Some(id) => {
                        let (mapped, unmapped) = counts
                            .get_mut(id)
                            .context("record reference ID out of range")?;
                        if record.flags().is_unmapped() {
                            *unmapped += 1;
                        } else {
                            *mapped += 1;
                        }
                    }
                    None => unplaced += 1,
                }
            }
            (counts, unplaced)
        }
    };
    for ((name, len), (mapped, unmapped)) in refs.iter().zip(counts) {
        writeln!(out, "{name}\t{len}\t{mapped}\t{unmapped}")?;
    }

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
