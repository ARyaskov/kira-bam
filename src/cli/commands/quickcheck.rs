use anyhow::Result;

use crate::cli::QuickcheckArgs;
use crate::quickcheck::check;

pub fn cmd_quickcheck(args: QuickcheckArgs) -> Result<()> {
    if let Err(failures) = check(&args.inputs, args.unrecognized_ok) {
        if args.verbose {
            for (path, error) in &failures.files {
                eprintln!("{}: {error:#}", path.display());
            }
        }
        return Err(failures.into());
    }
    Ok(())
}
