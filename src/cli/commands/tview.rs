use anyhow::Result;

use crate::cli::TviewArgs;
use crate::tview::run as run_tview;

pub fn cmd_tview(args: TviewArgs) -> Result<()> {
    run_tview(args)
}
