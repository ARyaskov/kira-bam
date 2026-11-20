use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::PathBuf;

use anyhow::Result;
use noodles_bam as bam;

use crate::cli::QuickcheckArgs;
use crate::io::BGZF_EOF;

#[derive(Debug)]
pub struct CheckFailures {
    pub files: Vec<(PathBuf, anyhow::Error)>,
}

impl std::fmt::Display for CheckFailures {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "quickcheck failed for {} file(s)", self.files.len())
    }
}

impl std::error::Error for CheckFailures {}

pub fn run(args: QuickcheckArgs) -> Result<()> {
    check(&args.inputs, args.unrecognized_ok)?;
    Ok(())
}

pub fn check(inputs: &[PathBuf], unrecognized_ok: bool) -> std::result::Result<(), CheckFailures> {
    let files = inputs
        .iter()
        .filter_map(|path| {
            check_one(path, unrecognized_ok)
                .err()
                .map(|err| (path.clone(), err))
        })
        .collect::<Vec<_>>();
    if files.is_empty() {
        Ok(())
    } else {
        Err(CheckFailures { files })
    }
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
