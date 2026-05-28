use anyhow::Result;

use crate::cli::FixmateArgs;
use crate::fixmate::run as run_fixmate;

pub fn cmd_fixmate(args: FixmateArgs) -> Result<()> {
    run_fixmate(args)
}
