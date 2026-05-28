use anyhow::Result;

use crate::cli::IndexArgs;
use crate::index::run as run_index;

pub fn cmd_index(args: IndexArgs) -> Result<()> {
    run_index(args)
}
