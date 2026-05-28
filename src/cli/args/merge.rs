use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct MergeArgs {
    /// Output BAM file.
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Sorted input BAM files (at least one).
    #[arg(value_name = "IN", required = true, num_args = 1..)]
    pub inputs: Vec<PathBuf>,

    /// Inputs are sorted by read name instead of coordinate.
    #[arg(short = 'n', long = "name-sort")]
    pub name_sort: bool,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Overwrite output if it exists.
    #[arg(short = 'f', long = "force")]
    pub force: bool,

    /// Write uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,

    /// Reference FASTA — required for CRAM input or output.
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,

    /// Do not append a @PG record to the output header.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}
