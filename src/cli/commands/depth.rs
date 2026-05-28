use anyhow::Result;

use crate::cli::DepthArgs;
use crate::depth::run as run_depth;

pub fn cmd_depth(args: DepthArgs) -> Result<()> {
    run_depth(args)
}
