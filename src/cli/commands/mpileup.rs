use anyhow::Result;

use crate::cli::MpileupArgs;
use crate::mpileup::run as run_mpileup;

pub fn cmd_mpileup(args: MpileupArgs) -> Result<()> {
    run_mpileup(args)
}
