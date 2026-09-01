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

    /// Number of temporary files to hash read names across when the input does
    /// not fit in memory. Every record of one name lands in the same file, so
    /// each can then be grouped independently. More files = less memory per
    /// file, at the cost of more open descriptors.
    #[arg(short = 'r', long = "buckets", default_value_t = 64)]
    pub buckets: usize,

    /// Approximate memory budget before spilling to `--tmpdir` (e.g. `768M`,
    /// `4G`, or `auto` for a fraction of total RAM).
    #[arg(short = 'm', long = "memory", default_value = "768M")]
    pub memory: String,

    /// Tmp dir.
    #[arg(short = 'T', long = "tmpdir")]
    pub tmpdir: Option<PathBuf>,

    /// Reference FASTA, required to decode CRAM input (`.fai` must exist
    /// alongside). Long-only: `-T` is the tmp dir here, as in samtools collate.
    #[arg(long = "reference")]
    pub reference: Option<PathBuf>,
}
