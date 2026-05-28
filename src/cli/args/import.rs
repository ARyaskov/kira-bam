use std::path::PathBuf;

use clap::Parser;

/// Convert FASTQ → unaligned BAM/SAM.
#[derive(Parser, Debug)]
pub struct ImportArgs {
    /// Input FASTQ for R1 (or single-end).
    #[arg(value_name = "R1")]
    pub r1: PathBuf,

    /// Input FASTQ for R2 (optional, paired-end).
    #[arg(value_name = "R2")]
    pub r2: Option<PathBuf>,

    /// Output BAM file.
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Sample name → @RG SM tag.
    #[arg(short = 's', long = "sample-name")]
    pub sample: Option<String>,

    /// Read group ID.
    #[arg(short = 'g', long = "rg-id", default_value = "1")]
    pub rg_id: String,

    /// Threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Write uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,

    /// Parse CASAVA `/1`/`/2` suffix and pull barcode (BC tag).
    #[arg(short = 'C', long = "casava")]
    pub casava: bool,

    /// Don't append @PG.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}
