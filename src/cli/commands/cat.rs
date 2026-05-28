use anyhow::Result;

use crate::cli::CatArgs;
use crate::cat::run as run_cat;

pub fn cmd_cat(args: CatArgs) -> Result<()> {
    run_cat(args)
}
