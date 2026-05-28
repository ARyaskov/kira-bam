use std::path::PathBuf;

use clap::Parser;

/// Add or replace @RG headers and per-record RG tags.
#[derive(Parser, Debug)]
pub struct AddReplaceRgArgs {
    /// Input BAM/SAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file.
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// New @RG line in `ID:foo\tSM:bar` form (one of -r or --rg-line is required).
    #[arg(short = 'r', long = "rg-line", required = true)]
    pub rg_line: String,

    /// Mode: orphan-only (default), overwrite_all, overwrite_orphans, lenient.
    #[arg(short = 'm', long = "mode", default_value = "orphan-only")]
    pub mode: String,

    /// Don't append a @PG record.
    #[arg(long = "no-PG")]
    pub no_pg: bool,

    /// Threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Write uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,
}
