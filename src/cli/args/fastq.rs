use std::path::PathBuf;

use clap::Parser;

/// Export BAM/SAM to FASTQ.
#[derive(Parser, Debug)]
pub struct FastqArgs {
    /// Input BAM/SAM file. `-` for stdin.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file for read 1 (or single-end when -2 not set).
    #[arg(short = '1')]
    pub read1: Option<PathBuf>,

    /// Output file for read 2.
    #[arg(short = '2')]
    pub read2: Option<PathBuf>,

    /// Output file for singletons (mate-missing reads).
    #[arg(short = 's', long = "singletons")]
    pub singletons: Option<PathBuf>,

    /// Output file for orphan reads (reads with FLAG_PAIRED but neither READ1 nor READ2).
    #[arg(short = '0', long = "orphans")]
    pub orphans: Option<PathBuf>,

    /// Output as FASTA instead of FASTQ.
    #[arg(long = "fasta")]
    pub fasta: bool,

    /// Include CASAVA-style `/1` `/2` suffix in read names.
    #[arg(short = 'N', long = "casava")]
    pub casava: bool,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Collate reads by name (assume input is not name-sorted).
    #[arg(long = "collate")]
    pub collate: bool,
}
