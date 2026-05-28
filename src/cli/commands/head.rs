use anyhow::Result;

use crate::cli::HeadArgs;
use crate::head::run as run_head;

pub fn cmd_head(args: HeadArgs) -> Result<()> {
    run_head(args)
}
