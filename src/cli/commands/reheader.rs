use anyhow::Result;

use crate::cli::ReheaderArgs;
use crate::reheader::run as run_reheader;

pub fn cmd_reheader(args: ReheaderArgs) -> Result<()> {
    run_reheader(args)
}
