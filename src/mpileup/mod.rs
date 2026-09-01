use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};

use crate::cli::MpileupArgs;
use crate::faidx::RefCache;
use crate::io::BamReader;
use crate::pileup::PileupIter;

pub fn run(args: MpileupArgs) -> Result<()> {
    if args.inputs.len() > 1 {
        eprintln!("[kira-bam mpileup] note: only the first BAM is processed in this version");
    }
    let input = &args.inputs[0];
    let reader =
        BamReader::open_with_reference(input, args.reference.as_deref()).context("open input")?;
    let header = reader.header().clone();

    // Indexed, one contig cached at a time — a whole-genome FASTA no longer
    // has to fit in RAM to call a pileup.
    let mut refs = match &args.reference {
        Some(p) => Some(RefCache::open(p)?),
        None => None,
    };
    // Contig names resolved once, not rebuilt as a String per column.
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| n.to_string())
        .collect();

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut pileup = PileupIter::new(reader, args.min_mapq, args.min_baseq);
    // One column buffer for the whole scan.
    let mut col: Vec<(u8, u8)> = Vec::new();
    while let Some((tid, pos)) = pileup.next_column_into(&mut col)? {
        if col.len() > args.max_depth as usize {
            col.truncate(args.max_depth as usize);
        }
        let name = ref_names.get(tid).map(String::as_str).unwrap_or("");
        let ref_base = match &mut refs {
            Some(cache) => cache
                .sequence(name)?
                .and_then(|s| s.get(pos.saturating_sub(1)))
                .copied()
                .unwrap_or(b'N'),
            None => b'N',
        };
        write!(out, "{name}\t{pos}\t{}\t{}", ref_base as char, col.len())?;
        // Read bases column with ref-base style markup (. for match, ACGT for mismatch).
        write!(out, "\t")?;
        for &(b, _q) in &col {
            let same = b.eq_ignore_ascii_case(&ref_base);
            if same {
                out.write_all(b".")?;
            } else {
                out.write_all(&[b.to_ascii_uppercase()])?;
            }
        }
        // Quality column.
        write!(out, "\t")?;
        for &(_b, q) in &col {
            // Saturating: the 0xFF "quality unavailable" sentinel overflows a
            // plain `+ 33`.
            out.write_all(&[q.saturating_add(33).min(126)])?;
        }
        writeln!(out)?;
    }
    Ok(())
}
