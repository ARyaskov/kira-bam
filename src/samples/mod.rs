use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;

use crate::cli::SamplesArgs;
use crate::io::BamReader;

pub fn run(args: SamplesArgs) -> Result<()> {
    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut seen: BTreeSet<String> = BTreeSet::new();
    for input in &args.inputs {
        let reader = BamReader::open(input).with_context(|| format!("open {}", input.display()))?;
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
    Ok(())
}
