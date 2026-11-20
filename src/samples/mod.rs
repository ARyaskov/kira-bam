use std::collections::BTreeSet;
use std::io::Write;

use anyhow::{Context, Result};
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;

use crate::cli::SamplesArgs;
use crate::io::BamReader;
use crate::io::BufferedOutput;

pub fn run(args: SamplesArgs) -> Result<()> {
    let mut out = BufferedOutput::create(args.output.as_deref())?;

    let mut seen: BTreeSet<String> = BTreeSet::new();
    for input in &args.inputs {
        let reader = BamReader::open_with_reference(input, args.reference.as_deref())
            .with_context(|| format!("open {}", input.display()))?;
        let header = reader.header();
        for (rg_id, rg_map) in header.read_groups() {
            let sm = rg_map
                .other_fields()
                .get(&rg_tag::SAMPLE)
                .map(|v| String::from_utf8_lossy(v).to_string())
                .unwrap_or_else(|| "-".to_string());
            if args.expand {
                writeln!(
                    out,
                    "{sm}\t{rg}\t{path}",
                    rg = rg_id,
                    path = input.display()
                )?;
            } else {
                seen.insert(sm);
            }
        }
    }
    if !args.expand {
        for s in seen {
            writeln!(out, "{s}")?;
        }
    }
    out.finish()
}
