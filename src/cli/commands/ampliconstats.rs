use anyhow::Result;

use crate::ampliconstats::run as run_ampliconstats;
use crate::cli::AmpliconstatsArgs;

pub fn cmd_ampliconstats(args: AmpliconstatsArgs) -> Result<()> {
    run_ampliconstats(args)
}
