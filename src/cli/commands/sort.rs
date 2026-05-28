use anyhow::Result;

use crate::cli::SortArgs;
use crate::sort::run as run_sort;

pub fn cmd_sort(args: SortArgs) -> Result<()> {
    run_sort(args)
}
