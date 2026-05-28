use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Write};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::op::Kind;
use rustc_hash::FxHasher;

use crate::cli::TviewArgs;

const RESET: &str = "\x1b[0m";
const BOLD: &str = "\x1b[1m";
const RED: &str = "\x1b[31m";
const GREEN: &str = "\x1b[32m";
const YELLOW: &str = "\x1b[33m";
const BLUE: &str = "\x1b[34m";
const GREY: &str = "\x1b[90m";

pub fn run(args: TviewArgs) -> Result<()> {
    let (chrom, start_1) = parse_position(&args.position).context("parse position")?;
    let end_1 = start_1 + (args.width as u64).saturating_sub(1);

    let region: Region = format!("{chrom}:{start_1}-{end_1}")
        .parse()
        .map_err(|e| anyhow::anyhow!("region parse: {e:?}"))?;

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&args.input)
        .context("open indexed BAM")?;
    let header = reader.read_header().context("read header")?;

    let query = reader.query(&header, &region)?;
    let mut reads: Vec<RecordBuf> = Vec::new();
    for r in query.records() {
        let r = r?;
        let buf = RecordBuf::try_from_alignment_record(&header, &r)?;
        reads.push(buf);
        if reads.len() >= args.depth.saturating_mul(4) {
            break;
        }
    }

    let refs = match &args.reference {
        Some(p) => load_fasta(p)?,
        None => HashMap::default(),
    };
    let ref_seq = refs.get(chrom.as_bytes());

    let mut stdout = std::io::stdout().lock();
    let color = !args.no_color;

    // Header line: position ruler
    writeln!(stdout, "{chrom}:{start_1}-{end_1}")?;
    for offset in 0..args.width {
        let pos = start_1 + offset as u64;
        if pos.is_multiple_of(10) && offset + 9 < args.width {
            write!(stdout, "|")?;
        } else if pos.is_multiple_of(5) {
            write!(stdout, ":")?;
        } else {
            write!(stdout, " ")?;
        }
    }
    writeln!(stdout)?;

    // Reference line
    if let Some(rs) = ref_seq {
        for offset in 0..args.width {
            let pos = start_1 + offset as u64;
            let idx = (pos.saturating_sub(1)) as usize;
            let b = if idx < rs.len() { rs[idx] } else { b'.' };
            if color {
                write!(stdout, "{BOLD}{}{}", b as char, RESET)?;
            } else {
                write!(stdout, "{}", b as char)?;
            }
        }
        writeln!(stdout)?;
    } else {
        writeln!(stdout, "{}", "N".repeat(args.width))?;
    }

    // Layout reads into stacked rows (greedy bin-pack)
    let mut rows: Vec<Vec<(u64, u64, Vec<u8>)>> = Vec::new(); // (start, end, bases)
    for rec in &reads {
        let aligned = render_read(rec, start_1, args.width);
        if aligned.is_none() {
            continue;
        }
        let (s, e, bases) = aligned.unwrap();
        let mut placed = false;
        for row in rows.iter_mut() {
            if row.last().map(|t| t.1 < s).unwrap_or(true) {
                row.push((s, e, bases.clone()));
                placed = true;
                break;
            }
        }
        if !placed && rows.len() < args.depth {
            rows.push(vec![(s, e, bases)]);
        }
    }

    for row in &rows {
        let mut line = vec![b' '; args.width];
        for (s, _e, bases) in row {
            let offset_in_row = s.saturating_sub(start_1) as usize;
            for (i, &b) in bases.iter().enumerate() {
                let idx = offset_in_row + i;
                if idx < args.width {
                    line[idx] = b;
                }
            }
        }
        if color {
            for &b in &line {
                let painted = paint(b, ref_seq, start_1, &line);
                stdout.write_all(painted.as_bytes())?;
            }
        } else {
            stdout.write_all(&line)?;
        }
        writeln!(stdout)?;
    }
    Ok(())
}

fn paint(b: u8, _ref_seq: Option<&Vec<u8>>, _start_1: u64, _line: &[u8]) -> String {
    let color = match b {
        b'A' | b'a' => GREEN,
        b'C' | b'c' => BLUE,
        b'G' | b'g' => YELLOW,
        b'T' | b't' => RED,
        b'N' | b'n' => GREY,
        b'-' => GREY,
        _ => "",
    };
    if color.is_empty() {
        (b as char).to_string()
    } else {
        format!("{color}{}{RESET}", b as char)
    }
}

fn render_read(rec: &RecordBuf, win_start_1: u64, win_width: usize) -> Option<(u64, u64, Vec<u8>)> {
    let start_1 = rec.alignment_start().map(|p| usize::from(p) as u64)?;
    let seq = rec.sequence().as_ref();
    let mut bases: Vec<u8> = Vec::new();
    let mut read_idx = 0usize;
    let mut ref_pos = start_1;
    let win_end_1 = win_start_1 + win_width as u64 - 1;
    for op in rec.cigar().as_ref().iter() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    if read_idx >= seq.len() {
                        break;
                    }
                    if ref_pos >= win_start_1 && ref_pos <= win_end_1 {
                        bases.push(seq[read_idx]);
                    }
                    ref_pos += 1;
                    read_idx += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                read_idx += len;
            }
            Kind::Deletion | Kind::Skip => {
                for _ in 0..len {
                    if ref_pos >= win_start_1 && ref_pos <= win_end_1 {
                        bases.push(b'-');
                    }
                    ref_pos += 1;
                }
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    if bases.is_empty() {
        return None;
    }
    let render_start = start_1.max(win_start_1);
    let render_end = render_start + bases.len() as u64 - 1;
    Some((render_start, render_end, bases))
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

fn parse_position(s: &str) -> Result<(String, u64)> {
    let (chr, pos_str) = s.split_once(':').context("position must be chr:pos")?;
    let pos: u64 = pos_str.parse().context("position must be an integer")?;
    Ok((chr.to_string(), pos))
}
