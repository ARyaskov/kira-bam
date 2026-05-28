use anyhow::Result;

use crate::addreplacerg::run as run_addreplacerg;
use crate::cli::AddReplaceRgArgs;

pub fn cmd_addreplacerg(args: AddReplaceRgArgs) -> Result<()> {
    run_addreplacerg(args)
}
