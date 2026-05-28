use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

use anyhow::{Context, Result};
use noodles_fasta as fasta;
use rustc_hash::FxHasher;

use crate::cli::MpileupArgs;
use crate::io::BamReader;
use crate::pileup::PileupIter;

pub fn run(args: MpileupArgs) -> Result<()> {
    if args.inputs.len() > 1 {
        eprintln!("[kira-bam mpileup] note: only the first BAM is processed in this version");
    }
    let input = &args.inputs[0];
    let reader = BamReader::open(input).context("open input")?;
    let header = reader.header().clone();

    let refs = match &args.reference {
        Some(p) => load_fasta(p)?,
        None => HashMap::default(),
    };

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut pileup = PileupIter::new(reader, args.min_mapq, args.min_baseq);
    while let Some((tid, pos, mut col)) = pileup.next_column()? {
        if col.len() > args.max_depth as usize {
            col.truncate(args.max_depth as usize);
        }
        let name: String = header
            .reference_sequences()
            .get_index(tid)
            .map(|(n, _)| n.to_string())
            .unwrap_or_default();
        let ref_base = refs
            .get(name.as_bytes())
            .and_then(|s| s.get(pos.saturating_sub(1)))
            .copied()
            .unwrap_or(b'N');
        write!(out, "{name}\t{pos}\t{}\t{}", ref_base as char, col.len())?;
        // Read bases column with ref-base style markup (. for match, ACGT for mismatch).
        write!(out, "\t")?;
        for &(b, _q) in &col {
            let same = b.to_ascii_uppercase() == ref_base.to_ascii_uppercase();
            if same {
                out.write_all(b".")?;
            } else {
                out.write_all(&[b.to_ascii_uppercase()])?;
            }
        }
        // Quality column.
        write!(out, "\t")?;
        for &(_b, q) in &col {
            out.write_all(&[(q + 33).min(126)])?;
        }
        writeln!(out)?;
    }
    Ok(())
}

fn load_fasta(
    path: &std::path::Path,
) -> Result<HashMap<Vec<u8>, Vec<u8>, std::hash::BuildHasherDefault<FxHasher>>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut reader = fasta::io::Reader::new(BufReader::with_capacity(1 << 20, f));
    let mut out: HashMap<Vec<u8>, Vec<u8>, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    for result in reader.records() {
        let r = result?;
        out.insert(r.name().to_vec(), r.sequence().as_ref().to_vec());
    }
    Ok(out)
}
