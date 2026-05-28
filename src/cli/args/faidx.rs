use std::path::PathBuf;

use clap::Parser;

/// Index FASTA (build .fai) and optionally extract regions.
#[derive(Parser, Debug)]
pub struct FaidxArgs {
    /// Reference FASTA.
    #[arg(value_name = "FASTA")]
    pub fasta: PathBuf,

    /// Region(s) to extract (`chr`, `chr:start-end`). If empty: just (re)build index.
    #[arg(value_name = "REGION")]
    pub regions: Vec<String>,

    /// Output file (default stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// FASTA line width when extracting.
    #[arg(long = "line-width", default_value_t = 60)]
    pub line_width: usize,
}
