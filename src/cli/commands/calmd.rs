use anyhow::Result;

use crate::cli::CalmdArgs;
use crate::calmd::run as run_calmd;

pub fn cmd_calmd(args: CalmdArgs) -> Result<()> {
    run_calmd(args)
}
