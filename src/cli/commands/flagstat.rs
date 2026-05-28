use anyhow::Result;

use crate::cli::FlagstatArgs;
use crate::flagstat::run as run_flagstat;

pub fn cmd_flagstat(args: FlagstatArgs) -> Result<()> {
    run_flagstat(args)
}
