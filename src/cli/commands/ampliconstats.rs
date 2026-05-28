use anyhow::Result;

use crate::cli::AmpliconstatsArgs;
use crate::ampliconstats::run as run_ampliconstats;

pub fn cmd_ampliconstats(args: AmpliconstatsArgs) -> Result<()> {
    run_ampliconstats(args)
}
