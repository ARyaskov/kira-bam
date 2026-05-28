use anyhow::Result;

use crate::cli::CoverageArgs;
use crate::coverage::run as run_coverage;

pub fn cmd_coverage(args: CoverageArgs) -> Result<()> {
    run_coverage(args)
}
