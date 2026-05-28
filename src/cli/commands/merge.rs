use anyhow::Result;

use crate::cli::MergeArgs;
use crate::merge::run as run_merge;

pub fn cmd_merge(args: MergeArgs) -> Result<()> {
    run_merge(args)
}
