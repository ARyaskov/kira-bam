use std::path::PathBuf;

use clap::Parser;

/// Print first N records of a BAM/SAM.
#[derive(Parser, Debug)]
pub struct HeadArgs {
    /// Input.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Number of records to print (default 5).
    #[arg(short = 'n', long = "num", default_value_t = 5)]
    pub n: u64,

    /// Number of header lines to print (samtools `-h`). 0 = all.
    #[arg(short = 'h', long = "header-lines", default_value_t = 0)]
    pub header_lines: usize,

    /// Output file (default stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,
}
