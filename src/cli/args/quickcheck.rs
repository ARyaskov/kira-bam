use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct QuickcheckArgs {
    /// Input BAM/SAM files.
    #[arg(value_name = "IN", required = true, num_args = 1..)]
    pub inputs: Vec<PathBuf>,

    /// Verbose: print one diagnostic per failing file.
    #[arg(short = 'v', long = "verbose")]
    pub verbose: bool,

    /// Treat unrecognized formats as OK (samtools default is failure).
    #[arg(short = 'u', long = "unrecognized-ok")]
    pub unrecognized_ok: bool,
}
