use anyhow::Result;

use crate::ampliconclip::run as run_ampliconclip;
use crate::cli::AmpliconclipArgs;

pub fn cmd_ampliconclip(args: AmpliconclipArgs) -> Result<()> {
    run_ampliconclip(args)
}
