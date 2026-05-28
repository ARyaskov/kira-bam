use anyhow::Result;

use crate::cli::ViewArgs;
use crate::view::run as run_view;

pub fn cmd_view(args: ViewArgs) -> Result<()> {
    run_view(args)
}
