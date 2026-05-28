use std::path::PathBuf;

use clap::Parser;

/// Comprehensive BAM statistics (mirrors `samtools stats`).
#[derive(Parser, Debug)]
pub struct StatsArgs {
    /// Input BAM/SAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file (default: stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Optional reference FASTA (enables GC and error-rate stats).
    #[arg(short = 'r', long = "reference")]
    pub reference: Option<PathBuf>,

    /// Restrict to regions in BED.
    #[arg(short = 't', long = "target-regions")]
    pub bed: Option<PathBuf>,

    /// Sample name annotation (samtools PR #1864).
    #[arg(long = "sample-name")]
    pub sample_name: Option<String>,

    /// Min MAPQ.
    #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,
}
