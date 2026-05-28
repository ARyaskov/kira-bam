use std::path::PathBuf;

use clap::Parser;

/// Coverage statistics by BED interval.
#[derive(Parser, Debug)]
pub struct BedcovArgs {
    /// BED file with intervals (3 cols: chr, start, end).
    #[arg(value_name = "BED")]
    pub bed: PathBuf,

    /// Input BAM(s). Each one becomes a separate column.
    #[arg(value_name = "IN", required = true, num_args = 1..)]
    pub inputs: Vec<PathBuf>,

    /// Minimum read MAPQ.
    #[arg(short = 'Q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Output (default stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,
}
