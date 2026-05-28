use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct CoverageArgs {
    /// BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Output file (default: stdout).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Region spec(s); empty = all references.
    #[arg(short = 'r', long = "region")]
    pub regions: Vec<String>,

    /// Minimum MAPQ.
    #[arg(short = 'q', long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Minimum base quality.
    #[arg(short = 'Q', long = "min-baseq", default_value_t = 0)]
    pub min_baseq: u8,

    /// Don't print header line.
    #[arg(short = 'H', long = "no-header")]
    pub no_header: bool,

    /// Comma-separated subset of columns to print (samtools issue #1664).
    /// Accepts: rname,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq
    #[arg(long = "columns")]
    pub columns: Option<String>,

    /// Number of threads.
    #[arg(short = '@', long = "threads", default_value_t = 4)]
    pub threads: usize,
}
