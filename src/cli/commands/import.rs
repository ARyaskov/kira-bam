use anyhow::Result;

use crate::cli::ImportArgs;
use crate::import::run as run_import;

pub fn cmd_import(args: ImportArgs) -> Result<()> {
    run_import(args)
}
