use std::path::PathBuf;

use clap::Parser;

/// Concatenate BAMs that share a header.
#[derive(Parser, Debug)]
pub struct CatArgs {
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
}
