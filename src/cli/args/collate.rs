use std::path::PathBuf;

use clap::Parser;

/// Group reads of the same qname together (approximate name-sort).
#[derive(Parser, Debug)]
pub struct CollateArgs {
    /// Input BAM/SAM.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output prefix (we'll write `<prefix>.bam`).
    #[arg(value_name = "OUT_PREFIX")]
    pub output: PathBuf,

    /// Fast mode: only collate primary alignments (skip secondary/supplementary).
    #[arg(short = 'f', long = "fast")]
    pub fast: bool,

    /// Threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Don't append a @PG record.
    #[arg(long = "no-PG")]
    pub no_pg: bool,

    /// Number of hash buckets (samtools `-r`). Larger = less spill but more RAM.
    #[arg(short = 'r', long = "buckets", default_value_t = 1_048_576)]
    pub buckets: usize,

    /// Tmp dir.
    #[arg(short = 'T', long = "tmpdir")]
    pub tmpdir: Option<PathBuf>,
}
