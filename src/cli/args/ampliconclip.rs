use std::path::PathBuf;

use clap::Parser;

/// Soft-clip reads at primer boundaries defined in a BED file.
#[derive(Parser, Debug)]
pub struct AmpliconclipArgs {
    /// Input BAM/SAM (coordinate-sorted).
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// BED file with primer intervals.
    #[arg(short = 'b', long = "bed", required = true)]
    pub bed: PathBuf,

    /// Output file.
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Hard-clip instead of soft-clip.
    #[arg(long = "hard-clip")]
    pub hard_clip: bool,

    /// Threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Write uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,

    /// Don't append a @PG record.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}
