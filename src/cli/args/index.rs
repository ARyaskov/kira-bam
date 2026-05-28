use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct IndexArgs {
    /// Input sorted BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output index path (default: <input>.bai or <input>.csi).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Force BAI format (default if depth ≤ 14 levels).
    #[arg(short = 'b', long = "bai")]
    pub bai: bool,

    /// Force CSI format (use for references with chromosomes > 2^29 bp).
    #[arg(short = 'c', long = "csi")]
    pub csi: bool,

    /// CSI min-shift parameter.
    #[arg(long = "min-shift", default_value_t = 14)]
    pub min_shift: u8,

    /// CSI depth parameter.
    #[arg(long = "depth", default_value_t = 5)]
    pub depth: u8,

    /// Number of threads for BGZF decompression during scan.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,
}
