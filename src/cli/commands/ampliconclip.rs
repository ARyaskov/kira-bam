use anyhow::Result;

use crate::cli::AmpliconclipArgs;
use crate::ampliconclip::run as run_ampliconclip;

pub fn cmd_ampliconclip(args: AmpliconclipArgs) -> Result<()> {
    run_ampliconclip(args)
}
