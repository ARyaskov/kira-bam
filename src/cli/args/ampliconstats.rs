use std::path::PathBuf;

use clap::Parser;

/// Compute statistics for amplicons defined in a BED file.
#[derive(Parser, Debug)]
pub struct AmpliconstatsArgs {
    /// Amplicon BED file.
    #[arg(value_name = "BED")]
    pub bed: PathBuf,

    /// Input BAM(s).
    #[arg(value_name = "IN", required = true, num_args = 1..)]
    pub inputs: Vec<PathBuf>,

    /// Output file.
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Min MAPQ.
    #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,
}
