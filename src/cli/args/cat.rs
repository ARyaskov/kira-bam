use std::path::PathBuf;

use clap::Parser;

/// Concatenate BAMs that share a header.
/// `-h FILE` takes the header from FILE here, as in samtools cat, so clap's
/// auto-generated short help flag is turned off and `--help` declared
/// explicitly.
#[derive(Parser, Debug)]
#[command(disable_help_flag = true)]
pub struct CatArgs {
    /// Print help.
    #[arg(long = "help", action = clap::ArgAction::Help)]
    pub help: Option<bool>,

    /// Output file.
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Inputs.
    #[arg(value_name = "IN", required = true, num_args = 1..)]
    pub inputs: Vec<PathBuf>,

    /// Compression level (-1 = default, 0 = uncompressed).
    #[arg(short = 'l', long = "compression-level", default_value_t = -1)]
    pub compression: i32,

    /// Use header from this file (overrides input #0's header).
    #[arg(short = 'h', long = "header-from")]
    pub header_from: Option<PathBuf>,

    /// Do not append a @PG record.
    #[arg(long = "no-PG")]
    pub no_pg: bool,

    /// Reference FASTA, required to decode CRAM input (`.fai` must exist alongside).
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,
}
