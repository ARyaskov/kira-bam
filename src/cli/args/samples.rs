use std::path::PathBuf;

use clap::Parser;

/// List sample identifiers (SM:) from @RG headers.
#[derive(Parser, Debug)]
pub struct SamplesArgs {
    /// Input BAM/SAM files.
    #[arg(value_name = "IN", required = true, num_args = 1..)]
    pub inputs: Vec<PathBuf>,

    /// Print one row per @RG (default: one row per sample).
    #[arg(short = 'X', long = "expand")]
    pub expand: bool,

    /// Output file.
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,
}
