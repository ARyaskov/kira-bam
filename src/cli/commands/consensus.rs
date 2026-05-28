use anyhow::Result;

use crate::cli::ConsensusArgs;
use crate::consensus::run as run_consensus;

pub fn cmd_consensus(args: ConsensusArgs) -> Result<()> {
    run_consensus(args)
}
