use std::path::PathBuf;

use clap::Parser;

/// Split a BAM by @RG or another scalar tag.
#[derive(Parser, Debug)]
pub struct SplitArgs {
    /// Input BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file format pattern; `%!` expands to the tag value (default `%*_%!.bam`).
    ///
    /// Missing tags use `unknown`; the literal string `unknown` uses `%75nknown`.
    /// Integers use the `%i` prefix, floats use `%f_` followed by hex bits, and hex strings use `%h`.
    #[arg(short = 'f', long = "format", default_value = "%*_%!.bam")]
    pub format: String,

    /// Scalar tag to split on (default `RG`).
    #[arg(short = 'd', long = "tag", default_value = "RG")]
    pub tag: String,

    /// Don't append a @PG record.
    #[arg(long = "no-PG")]
    pub no_pg: bool,

    /// Accepted for CLI compatibility. This command runs single-threaded.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Cap on simultaneously open output files. Every distinct tag value needs
    /// its own writer (an fd plus a ~1 MiB buffer), so an unbounded split on a
    /// high-cardinality tag such as `CB` exhausts the process fd limit.
    #[arg(long = "max-outputs", default_value_t = 256)]
    pub max_outputs: usize,

    /// Reference FASTA, required to decode CRAM input (`.fai` must exist alongside).
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,
}
