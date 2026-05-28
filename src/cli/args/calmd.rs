use std::path::PathBuf;

use clap::Parser;

/// Recompute MD/NM tags by comparing reads to a reference FASTA.
#[derive(Parser, Debug)]
pub struct CalmdArgs {
    /// Input BAM/SAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Reference FASTA (must be indexed with .fai).
    #[arg(value_name = "REF")]
    pub reference: PathBuf,

    /// Output file (default: stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Output BAM.
    #[arg(short = 'b', long = "bam")]
    pub bam: bool,

    /// Write uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,

    /// Adjust BQ (BAQ-style) — currently a no-op stub.
    #[arg(short = 'e', long = "extended-baq")]
    pub extended_baq: bool,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Do not append a @PG record.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}
