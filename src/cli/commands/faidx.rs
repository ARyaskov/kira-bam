use anyhow::Result;

use crate::cli::FaidxArgs;
use crate::faidx::run as run_faidx;

pub fn cmd_faidx(args: FaidxArgs) -> Result<()> {
    run_faidx(args)
}
