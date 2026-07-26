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

    /// Optical-duplicate pixel distance (samtools `-d`). When > 0, kira-bam
    /// parses Illumina cluster coordinates (tile:x:y) from read names and reports
    /// the optical subset of PCR duplicates in the `-s` stats file. 0 disables it.
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

    /// Also flag the supplementary and secondary alignments of a duplicate
    /// template (samtools `-S`). Off by default, as in samtools; it costs one
    /// read-name entry per duplicate template. Without it those records keep a
    /// clear 0x400 while their primary is marked, so `-F 0x400` still lets them
    /// through.
    #[arg(short = 'S', long = "supplementary")]
    pub supplementary: bool,

    /// Do not append a @PG record to the output header.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}
