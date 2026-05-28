use anyhow::Result;

use crate::cli::AddReplaceRgArgs;
use crate::addreplacerg::run as run_addreplacerg;

pub fn cmd_addreplacerg(args: AddReplaceRgArgs) -> Result<()> {
    run_addreplacerg(args)
}
