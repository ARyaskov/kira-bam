use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct FlagstatArgs {
    /// Input BAM/SAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Number of threads (BGZF decompression).
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Output format: `text` (default) or `json`.
    #[arg(short = 'O', long = "output-fmt", default_value = "text")]
    pub output_fmt: String,
}
