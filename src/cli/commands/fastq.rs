use anyhow::Result;

use crate::cli::FastqArgs;
use crate::fastq::run as run_fastq;

pub fn cmd_fastq(args: FastqArgs) -> Result<()> {
    run_fastq(args)
}
