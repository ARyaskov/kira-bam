use anyhow::Result;

use crate::cli::SplitArgs;
use crate::split::run as run_split;

pub fn cmd_split(args: SplitArgs) -> Result<()> {
    run_split(args)
}
