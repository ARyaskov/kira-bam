use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct IdxstatsArgs {
    /// Indexed BAM file (.bai or .csi expected alongside).
    #[arg(value_name = "IN")]
    pub input: PathBuf,
}
