use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

use anyhow::Result;
use noodles_bam as bam;

use crate::cli::QuickcheckArgs;

/// BGZF EOF marker — empty BGZF block (28 bytes).
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

pub fn run(args: QuickcheckArgs) -> Result<()> {
    let mut bad_files: Vec<(std::path::PathBuf, String)> = Vec::new();
    for p in &args.inputs {
        match check_one(p, args.unrecognized_ok) {
            Ok(()) => {}
            Err(e) => bad_files.push((p.clone(), format!("{e:#}"))),
        }
    }
    if args.verbose {
        for (p, e) in &bad_files {
            eprintln!("{}: {e}", p.display());
        }
    }
    // samtools convention: exit non-zero if any bad files.
    if !bad_files.is_empty() {
        std::process::exit(1);
    }
    Ok(())
}

fn check_one(path: &std::path::Path, unrecognized_ok: bool) -> Result<()> {
    let mut f = File::open(path)?;
    let mut magic = [0u8; 4];
    let n = f.read(&mut magic)?;
    if n < 2 {
        anyhow::bail!("file too short");
    }
    let is_bgzf = magic[0] == 0x1f && magic[1] == 0x8b;
    if !is_bgzf {
        // Could be SAM. samtools quickcheck only validates BAM/CRAM.
        if unrecognized_ok {
            return Ok(());
        }
        anyhow::bail!("not a BGZF file (no gzip magic)");
    }

    // Validate header parses.
    let file2 = File::open(path)?;
    let mut reader = bam::io::Reader::new(file2);
    let _header = reader.read_header()?;

    // EOF marker check.
    let total_len = f.metadata()?.len();
    if total_len < BGZF_EOF.len() as u64 {
        anyhow::bail!("file shorter than EOF marker");
    }
    let mut tail = [0u8; 28];
    f.seek(SeekFrom::End(-28))?;
    f.read_exact(&mut tail)?;
    if tail != BGZF_EOF {
        anyhow::bail!("missing BGZF EOF marker — likely truncated");
    }
    Ok(())
}
