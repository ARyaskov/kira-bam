use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct DepthArgs {
    /// Indexed BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file (default: stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Output every position, including zero-coverage (samtools -a).
    #[arg(short = 'a', long = "all-positions")]
    pub all_positions: bool,

    /// BED file restricting positions reported.
    #[arg(short = 'b', long = "bed")]
    pub bed: Option<PathBuf>,

    /// Region spec, repeatable.
    #[arg(short = 'r', long = "region")]
    pub regions: Vec<String>,

    /// Minimum read MAPQ.
    #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Minimum base quality.
    #[arg(short = 'Q', long = "min-baseq", default_value_t = 0)]
    pub min_baseq: u8,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,
}
