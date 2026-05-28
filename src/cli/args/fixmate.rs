use std::path::PathBuf;

use clap::Parser;

/// Fix mate-related fields after name-sort: mate position, mate ref, MC/MS/MQ tags,
/// and 0x2 (proper pair) bit. Input must be queryname-sorted.
#[derive(Parser, Debug)]
pub struct FixmateArgs {
    /// Name-sorted input.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file.
    #[arg(value_name = "OUT")]
    pub output: PathBuf,

    /// Remove unmapped reads and secondary/supplementary alignments (samtools -r).
    #[arg(short = 'r', long = "remove-unpaired")]
    pub remove_unpaired: bool,

    /// Add MC (mate CIGAR) and MS (mate score) tags (samtools -m).
    #[arg(short = 'm', long = "add-mc-ms")]
    pub add_mc_ms: bool,

    /// Recompute CIGAR alignment from records (samtools -c).
    #[arg(short = 'c', long = "recompute-cigar")]
    pub recompute_cigar: bool,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Write uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,

    /// Do not append a @PG record to the output header.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}
