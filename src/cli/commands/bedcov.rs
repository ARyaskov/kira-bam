use anyhow::Result;

use crate::bedcov::run as run_bedcov;
use crate::cli::BedcovArgs;

pub fn cmd_bedcov(args: BedcovArgs) -> Result<()> {
    run_bedcov(args)
}
