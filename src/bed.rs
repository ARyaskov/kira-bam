use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Result, ensure};

#[derive(Debug, Clone)]
pub(crate) struct Interval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

pub(crate) fn read(path: &std::path::Path) -> Result<Vec<Interval>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    for (index, line) in r.lines().enumerate() {
        let location = || format!("{}: line {}", path.display(), index + 1);
        let line = line.with_context(location)?;
        if line.trim().is_empty()
            || line.starts_with('#')
            || line.starts_with("track ")
            || line.starts_with("browser ")
        {
            continue;
        }
        out.push(parse_interval(&line).with_context(location)?);
    }
    Ok(out)
}

fn parse_interval(line: &str) -> Result<Interval> {
    let mut fields = line.split('\t');
    let chrom = fields
        .next()
        .filter(|field| !field.is_empty())
        .context("missing chromosome")?;
    let start = fields
        .next()
        .context("missing BED start")?
        .parse()
        .context("invalid BED start")?;
    let end = fields
        .next()
        .context("missing BED end")?
        .parse()
        .context("invalid BED end")?;
    ensure!(start <= end, "BED start exceeds end");
    ensure!(
        usize::try_from(end).is_ok(),
        "BED end exceeds addressable coordinates"
    );
    Ok(Interval {
        chrom: chrom.to_owned(),
        start,
        end,
    })
}
