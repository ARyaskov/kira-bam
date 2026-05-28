use anyhow::Result;

use crate::cli::MarkdupArgs;
use crate::markdup::run as run_markdup;

pub fn cmd_markdup(args: MarkdupArgs) -> Result<()> {
    run_markdup(args)
}
