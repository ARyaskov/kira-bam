use anyhow::Result;

use crate::cli::StatsArgs;
use crate::stats::run as run_stats;

pub fn cmd_stats(args: StatsArgs) -> Result<()> {
    run_stats(args)
}
