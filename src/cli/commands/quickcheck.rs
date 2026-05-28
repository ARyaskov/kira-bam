use anyhow::Result;

use crate::cli::QuickcheckArgs;
use crate::quickcheck::run as run_quickcheck;

pub fn cmd_quickcheck(args: QuickcheckArgs) -> Result<()> {
    run_quickcheck(args)
}
