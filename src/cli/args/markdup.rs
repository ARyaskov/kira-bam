use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct MarkdupArgs {
    /// Sorted input BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output BAM file.
    #[arg(value_name = "OUT")]
    pub output: PathBuf,

    /// Remove duplicates instead of just marking the 0x400 flag.
    #[arg(short = 'r', long = "remove-dups")]
    pub remove: bool,

    /// Treat optical duplicates separately (samtools `-d` distance).
    #[arg(short = 'd', long = "optical-distance", default_value_t = 0)]
    pub optical_distance: u32,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Statistics output file.
    #[arg(short = 's', long = "stats")]
    pub stats: Option<PathBuf>,

    /// Write uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,

    /// Reference FASTA — required when input or output is CRAM.
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,

    /// Barcode tag (e.g. `RX`) used to group duplicates by UMI / molecular index.
    #[arg(long = "barcode-tag")]
    pub barcode_tag: Option<String>,

    /// Ancient-DNA mode: collapse forward and reverse strands together.
    #[arg(long = "mode-ancient")]
    pub mode_ancient: bool,

    /// Do not append a @PG record to the output header.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}
