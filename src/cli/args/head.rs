use std::path::PathBuf;

use clap::Parser;

/// Print first N records of a BAM/SAM.
/// `-h` is the header-line count here, as in samtools head, so clap's
/// auto-generated short help flag is turned off and `--help` declared
/// explicitly.
#[derive(Parser, Debug)]
#[command(disable_help_flag = true)]
pub struct HeadArgs {
    /// Print help.
    #[arg(long = "help", action = clap::ArgAction::Help)]
    pub help: Option<bool>,

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

    /// Reference FASTA, required to decode CRAM input (`.fai` must exist alongside).
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,
}
