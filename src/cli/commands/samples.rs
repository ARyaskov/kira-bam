use anyhow::Result;

use crate::cli::SamplesArgs;
use crate::samples::run as run_samples;

pub fn cmd_samples(args: SamplesArgs) -> Result<()> {
    run_samples(args)
}
