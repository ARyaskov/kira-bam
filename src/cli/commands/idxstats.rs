use anyhow::Result;

use crate::cli::IdxstatsArgs;
use crate::idxstats::run as run_idxstats;

pub fn cmd_idxstats(args: IdxstatsArgs) -> Result<()> {
    run_idxstats(args)
}
