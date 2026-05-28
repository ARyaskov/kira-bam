use anyhow::Result;

use crate::cli::FqidxArgs;
use crate::fqidx::run as run_fqidx;

pub fn cmd_fqidx(args: FqidxArgs) -> Result<()> {
    run_fqidx(args)
}
