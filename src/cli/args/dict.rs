use std::path::PathBuf;

use clap::Parser;

/// Build a sequence dictionary (.dict) from a FASTA.
#[derive(Parser, Debug)]
pub struct DictArgs {
    /// Reference FASTA.
    #[arg(value_name = "FASTA")]
    pub fasta: PathBuf,

    /// Output (default: <fasta>.dict on stdout if `-`).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Assembly identifier (AS field).
    #[arg(short = 'a', long = "assembly")]
    pub assembly: Option<String>,

    /// Species (SP field).
    #[arg(short = 's', long = "species")]
    pub species: Option<String>,

    /// URI field (UR).
    #[arg(short = 'u', long = "uri")]
    pub uri: Option<String>,
}
