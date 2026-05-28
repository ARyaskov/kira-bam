use std::path::PathBuf;

use clap::Parser;

/// Replace the header of a BAM file without rewriting the body.
#[derive(Parser, Debug)]
pub struct ReheaderArgs {
    /// New header (SAM format).
    #[arg(value_name = "HEADER")]
    pub header: PathBuf,

    /// Input BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output BAM file. Required (unlike samtools, we never overwrite in place).
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Don't append a @PG record.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}
