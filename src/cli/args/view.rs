use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct ViewArgs {
    /// Input file (SAM/BAM, or `-` for stdin).
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Optional region filters, `chr:start-end` (requires indexed BAM).
    #[arg(value_name = "REGION")]
    pub regions: Vec<String>,

    /// Region(s) to EXCLUDE (samtools open PR #669). Repeatable.
    #[arg(short = 'v', long = "exclude-region")]
    pub exclude_regions: Vec<String>,

    /// Output to FILE (default: stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Output BAM (default for non-stdout when input is BAM).
    #[arg(short = 'b', long = "bam")]
    pub bam: bool,

    /// Output SAM (default when writing to stdout).
    #[arg(short = 'S', long = "sam")]
    pub sam: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,

    /// Output CRAM (requires --reference).
    #[arg(short = 'C', long = "cram")]
    pub cram: bool,

    /// Reference FASTA (required for CRAM I/O). `.fai` must exist alongside.
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,

    /// Include header in SAM output.
    #[arg(short = 'h', long = "with-header")]
    pub with_header: bool,

    /// Only output header.
    #[arg(short = 'H', long = "header-only")]
    pub header_only: bool,

    /// Require all of the FLAGS in INT to be present (accepts decimal or 0x-hex).
    #[arg(short = 'f', long = "require-flags", default_value_t = 0, value_parser = parse_flags)]
    pub require_flags: u16,

    /// Filter out reads with any of the FLAGS in INT (accepts decimal or 0x-hex).
    #[arg(short = 'F', long = "filter-flags", default_value_t = 0, value_parser = parse_flags)]
    pub filter_flags: u16,

    /// Minimum MAPQ.
    #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,

    /// Count alignments instead of printing them.
    #[arg(short = 'c', long = "count")]
    pub count: bool,

    /// Drop all aux tags from output (samtools view -X).
    #[arg(short = 'X', long = "drop-tags")]
    pub drop_tags: bool,

    /// Do not append a @PG record to the output header.
    #[arg(long = "no-PG")]
    pub no_pg: bool,
}

fn parse_flags(s: &str) -> Result<u16, String> {
    let s = s.trim();
    if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u16::from_str_radix(hex, 16).map_err(|e| e.to_string())
    } else {
        s.parse::<u16>().map_err(|e| e.to_string())
    }
}
