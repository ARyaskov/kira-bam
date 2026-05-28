use anyhow::Result;

use crate::cli::CollateArgs;
use crate::collate::run as run_collate;

pub fn cmd_collate(args: CollateArgs) -> Result<()> {
    run_collate(args)
}
