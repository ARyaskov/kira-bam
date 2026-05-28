use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};

use crate::cli::ConsensusArgs;
use crate::io::BamReader;
use crate::pileup::PileupIter;

pub fn run(args: ConsensusArgs) -> Result<()> {
    let reader = BamReader::open(&args.input).context("open input")?;
    let header = reader.header().clone();

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut pileup = PileupIter::new(reader, args.min_mapq, args.min_baseq);

    let mut current_tid: Option<usize> = None;
    let mut last_pos: usize = 0;
    let mut line_pos: usize = 0;

    while let Some((tid, pos, col)) = pileup.next_column()? {
        if current_tid != Some(tid) {
            if current_tid.is_some() {
                writeln!(out)?;
            }
            let name = header
                .reference_sequences()
                .get_index(tid)
                .map(|(n, _)| n.to_string())
                .unwrap_or_default();
            writeln!(out, ">{name}")?;
            current_tid = Some(tid);
            last_pos = 0;
            line_pos = 0;
        }
        // Fill gaps with N.
        while last_pos + 1 < pos {
            write_base(&mut out, b'N', &mut line_pos, args.line_width)?;
            last_pos += 1;
        }
        let mut depth = col.len() as u32;
        if depth > args.max_depth {
            depth = args.max_depth;
        }
        let base = majority(&col, args.call_fraction);
        write_base(&mut out, base, &mut line_pos, args.line_width)?;
        last_pos = pos;
        let _ = depth;
    }
    if current_tid.is_some() {
        writeln!(out)?;
    }
    Ok(())
}

fn write_base(
    out: &mut Box<dyn Write>,
    b: u8,
    line_pos: &mut usize,
    width: usize,
) -> Result<()> {
    out.write_all(&[b])?;
    *line_pos += 1;
    if width > 0 && *line_pos >= width {
        out.write_all(b"\n")?;
        *line_pos = 0;
    }
    Ok(())
}

fn majority(col: &[(u8, u8)], threshold: f32) -> u8 {
    let mut counts = [0u32; 6]; // A C G T N other
    for &(b, _) in col {
        let i = match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 4,
            _ => 5,
        };
        counts[i] += 1;
    }
    let total: u32 = counts.iter().sum();
    if total == 0 {
        return b'N';
    }
    let (best_i, best) = counts
        .iter()
        .enumerate()
        .max_by_key(|&(_, &c)| c)
        .unwrap_or((4, &0));
    if (*best as f32) / (total as f32) < threshold {
        return b'N';
    }
    match best_i {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}
