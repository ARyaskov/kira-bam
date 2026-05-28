use std::path::PathBuf;

use clap::Parser;

/// Index FASTQ and optionally extract reads by name.
#[derive(Parser, Debug)]
pub struct FqidxArgs {
    /// FASTQ file.
    #[arg(value_name = "FASTQ")]
    pub fastq: PathBuf,

    /// Read names to extract (each becomes a 4-line FASTQ record on output).
    #[arg(value_name = "NAME")]
    pub names: Vec<String>,

    /// Output (default stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,
}
