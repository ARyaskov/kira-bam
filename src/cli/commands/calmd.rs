use anyhow::Result;

use crate::calmd::run as run_calmd;
use crate::cli::CalmdArgs;

pub fn cmd_calmd(args: CalmdArgs) -> Result<()> {
    run_calmd(args)
}
