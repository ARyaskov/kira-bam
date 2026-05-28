use std::path::PathBuf;

use clap::Parser;

/// Split a BAM by @RG (or any tag value).
#[derive(Parser, Debug)]
pub struct SplitArgs {
    /// Input BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file format pattern; `%!` expands to the tag value (default `%*_%!.bam`).
    #[arg(short = 'f', long = "format", default_value = "%*_%!.bam")]
    pub format: String,

    /// Tag to split on (default `RG`).
    #[arg(short = 'd', long = "tag", default_value = "RG")]
    pub tag: String,

    /// Don't append a @PG record.
    #[arg(long = "no-PG")]
    pub no_pg: bool,

    /// Threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,
}
