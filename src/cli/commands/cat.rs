use anyhow::Result;

use crate::cat::run as run_cat;
use crate::cli::CatArgs;

pub fn cmd_cat(args: CatArgs) -> Result<()> {
    run_cat(args)
}
