use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct MpileupArgs {
    /// Input BAM (sorted).
    #[arg(value_name = "IN", required = true)]
    pub inputs: Vec<PathBuf>,

    /// Output file.
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Optional reference FASTA.
    #[arg(short = 'f', long = "fasta-ref")]
    pub reference: Option<PathBuf>,

    /// Region(s) to scan.
    #[arg(short = 'r', long = "region")]
    pub regions: Vec<String>,

    /// Minimum MAPQ.
    #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Minimum base quality.
    #[arg(short = 'Q', long = "min-baseq", default_value_t = 13)]
    pub min_baseq: u8,

    /// Maximum depth per position.
    #[arg(short = 'd', long = "max-depth", default_value_t = 8000)]
    pub max_depth: u32,

    /// Output a row for every position even if zero coverage.
    #[arg(short = 'a', long = "all-positions")]
    pub all_positions: bool,
}
