use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct DepthArgs {
    /// Indexed BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file (default: stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Output every position of every referenced contig, including
    /// zero-coverage ones (samtools `-a`). Repeat (`-aa`) to also emit contigs
    /// with no reads at all.
    #[arg(short = 'a', long = "all-positions", action = clap::ArgAction::Count)]
    pub all_positions: u8,

    /// Output absolutely all positions, including unused reference sequences
    /// (samtools `--aa`; same as `-aa`).
    #[arg(long = "aa")]
    pub aa: bool,

    /// BED file restricting positions reported.
    #[arg(short = 'b', long = "bed")]
    pub bed: Option<PathBuf>,

    /// Region spec, repeatable.
    #[arg(short = 'r', long = "region")]
    pub regions: Vec<String>,

    /// Minimum read MAPQ.
    #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Filter out reads with any of these FLAGS (decimal or 0x-hex). The
    /// default is samtools' — UNMAP|SECONDARY|QCFAIL|DUP — without which
    /// duplicates and secondary alignments inflate every depth reported.
    #[arg(short = 'G', long = "excl-flags", default_value = "0x704", value_parser = parse_flags)]
    pub excl_flags: u16,

    /// Require all of these FLAGS to be present (decimal or 0x-hex).
    #[arg(short = 'g', long = "incl-flags", default_value_t = 0, value_parser = parse_flags)]
    pub incl_flags: u16,

    /// Minimum base quality.
    #[arg(short = 'Q', long = "min-baseq", default_value_t = 0)]
    pub min_baseq: u8,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Reference FASTA, required to decode CRAM input (`.fai` must exist alongside).
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,
}

fn parse_flags(s: &str) -> Result<u16, String> {
    let s = s.trim();
    if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u16::from_str_radix(hex, 16).map_err(|e| e.to_string())
    } else {
        s.parse::<u16>().map_err(|e| e.to_string())
    }
}
