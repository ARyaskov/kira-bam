use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct SortArgs {
    /// Input BAM/SAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output BAM file (default: stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Sort by read name (queryname) instead of coordinate.
    #[arg(short = 'n', long = "name-sort")]
    pub name_sort: bool,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 8)]
    pub threads: usize,

    /// Approximate per-thread memory budget for in-memory chunks (e.g. `768M`, `4G`).
    #[arg(short = 'm', long = "memory", default_value = "768M")]
    pub memory: String,

    /// Directory for temporary files (default: system tmp).
    #[arg(short = 'T', long = "tmpdir")]
    pub tmpdir: Option<PathBuf>,

    /// Write uncompressed BAM (level 0).
    #[arg(short = 'u', long = "uncompressed")]
    pub uncompressed: bool,

    /// Reference FASTA — only needed when input or output is CRAM.
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,

    /// Drop records without all of these flag bits set (samtools issue #912).
    #[arg(short = 'f', long = "require-flags", default_value_t = 0, value_parser = parse_flags)]
    pub require_flags: u16,

    /// Drop records with any of these flag bits set.
    #[arg(short = 'F', long = "filter-flags", default_value_t = 0, value_parser = parse_flags)]
    pub filter_flags: u16,

    /// Do not append a @PG record to the output header.
    #[arg(long = "no-PG")]
    pub no_pg: bool,

    /// Mark PCR/optical duplicates during the sort pass (FLAG 0x400).
    ///
    /// Equivalent to piping `kira-bam sort … | kira-bam markdup -`, but
    /// performed in-memory after the sort step — avoids one full BAM
    /// read+write roundtrip. Only effective when the input fits in a
    /// single in-memory chunk (i.e. `-m` is large enough); spilled sorts
    /// fall back to the standalone two-pass markdup.
    #[arg(long = "markdup")]
    pub markdup: bool,

    /// SAM tag (2 chars) holding cell barcode for single-cell-aware
    /// dedup. Only used with `--markdup`.
    #[arg(long = "markdup-barcode-tag", value_name = "TAG")]
    pub markdup_barcode_tag: Option<String>,

    /// "Ancient DNA" mode — ignore strand in dup key. Only with `--markdup`.
    #[arg(long = "markdup-mode-ancient")]
    pub markdup_mode_ancient: bool,
}

fn parse_flags(s: &str) -> Result<u16, String> {
    let s = s.trim();
    if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u16::from_str_radix(hex, 16).map_err(|e| e.to_string())
    } else {
        s.parse::<u16>().map_err(|e| e.to_string())
    }
}
