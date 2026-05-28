use anyhow::Result;

use crate::cli::BedcovArgs;
use crate::bedcov::run as run_bedcov;

pub fn cmd_bedcov(args: BedcovArgs) -> Result<()> {
    run_bedcov(args)
}
