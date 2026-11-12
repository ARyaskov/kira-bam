use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct ConsensusArgs {
    /// Input BAM (must be coordinate-sorted).
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output FASTA file.
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Minimum read MAPQ.
    #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Minimum base quality.
    #[arg(short = 'Q', long = "min-baseq", default_value_t = 13)]
    pub min_baseq: u8,

    /// Minimum observed depth to call a base (samtools consensus `-d`).
    /// Positions with fewer observations are emitted as `N`. 0 disables the
    /// minimum.
    #[arg(short = 'd', long = "min-depth", default_value_t = 1)]
    pub min_depth: u32,

    /// Threshold fraction for majority call (0.5–1.0). Below this → N.
    #[arg(long = "call-fraction", default_value_t = 0.5)]
    pub call_fraction: f32,

    /// Line wrap width in FASTA output.
    #[arg(long = "line-width", default_value_t = 60)]
    pub line_width: usize,

    /// Reference FASTA, required to decode CRAM input (`.fai` must exist alongside).
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,
}
