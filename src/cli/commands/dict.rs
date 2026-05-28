use anyhow::Result;

use crate::cli::DictArgs;
use crate::dict::run as run_dict;

pub fn cmd_dict(args: DictArgs) -> Result<()> {
    run_dict(args)
}
