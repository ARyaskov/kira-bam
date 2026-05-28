use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

use anyhow::{Context, Result};

use crate::cli::FaidxArgs;

#[derive(Debug, Clone)]
struct FaiEntry {
    name: String,
    length: u64,
    offset: u64,
    line_blen: u64,
    line_len: u64,
}

pub fn run(args: FaidxArgs) -> Result<()> {
    let fai_path = build_fai_path(&args.fasta);
    let fai = if !fai_path.exists() {
        let entries = build_fai(&args.fasta).context("build .fai")?;
        write_fai(&fai_path, &entries).context("write .fai")?;
        entries
    } else {
        read_fai(&fai_path).context("read .fai")?
    };

    if args.regions.is_empty() {
        return Ok(());
    }

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    let mut fasta_f = File::open(&args.fasta)?;

    for region in &args.regions {
        let (name, start, end) = parse_region(region);
        let entry = match fai.iter().find(|e| e.name == name) {
            Some(e) => e,
            None => {
                eprintln!("[kira-bam faidx] reference {name} not found");
                continue;
            }
        };
        let start = start.unwrap_or(1);
        let end = end.unwrap_or(entry.length);
        writeln!(out, ">{name}:{start}-{end}")?;
        let mut remaining = (end - start + 1) as usize;
        let mut pos_in_seq = (start - 1) as u64;
        while remaining > 0 {
            let line_idx = pos_in_seq / entry.line_blen;
            let in_line = pos_in_seq % entry.line_blen;
            let file_offset = entry.offset + line_idx * entry.line_len + in_line;
            let chunk = (entry.line_blen - in_line).min(remaining as u64) as usize;
            fasta_f.seek(SeekFrom::Start(file_offset))?;
            let mut buf = vec![0u8; chunk];
            fasta_f.read_exact(&mut buf)?;
            for chunk_part in buf.chunks(args.line_width.max(1)) {
                out.write_all(chunk_part)?;
                out.write_all(b"\n")?;
            }
            pos_in_seq += chunk as u64;
            remaining -= chunk;
        }
    }
    Ok(())
}

fn build_fai_path(fasta: &std::path::Path) -> std::path::PathBuf {
    let mut s = fasta.as_os_str().to_os_string();
    s.push(".fai");
    std::path::PathBuf::from(s)
}

fn build_fai(fasta: &std::path::Path) -> Result<Vec<FaiEntry>> {
    let f = File::open(fasta).context("open FASTA")?;
    let mut reader = BufReader::new(f);

    let mut entries: Vec<FaiEntry> = Vec::new();
    let mut current: Option<FaiEntry> = None;
    let mut file_pos: u64 = 0;
    let mut buf = String::new();

    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            break;
        }
        let trimmed_no_newline = buf.trim_end_matches('\n').trim_end_matches('\r');
        let line_bytes = n as u64;
        let line_blen = trimmed_no_newline.len() as u64;

        if buf.starts_with('>') {
            if let Some(e) = current.take() {
                entries.push(e);
            }
            let name = trimmed_no_newline[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            file_pos += line_bytes;
            current = Some(FaiEntry {
                name,
                length: 0,
                offset: file_pos,
                line_blen: 0,
                line_len: 0,
            });
        } else if let Some(c) = current.as_mut() {
            if c.line_blen == 0 {
                c.line_blen = line_blen;
                c.line_len = line_bytes;
            }
            c.length += line_blen;
            file_pos += line_bytes;
        } else {
            file_pos += line_bytes;
        }
    }
    if let Some(e) = current.take() {
        entries.push(e);
    }
    Ok(entries)
}

fn write_fai(path: &std::path::Path, entries: &[FaiEntry]) -> Result<()> {
    let mut f = BufWriter::new(File::create(path)?);
    for e in entries {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            e.name, e.length, e.offset, e.line_blen, e.line_len
        )?;
    }
    Ok(())
}

fn read_fai(path: &std::path::Path) -> Result<Vec<FaiEntry>> {
    let f = File::open(path)?;
    let mut entries = Vec::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }
        entries.push(FaiEntry {
            name: parts[0].to_string(),
            length: parts[1].parse().unwrap_or(0),
            offset: parts[2].parse().unwrap_or(0),
            line_blen: parts[3].parse().unwrap_or(0),
            line_len: parts[4].parse().unwrap_or(0),
        });
    }
    Ok(entries)
}

fn parse_region(s: &str) -> (String, Option<u64>, Option<u64>) {
    if let Some(colon) = s.find(':') {
        let (name, rest) = s.split_at(colon);
        let rest = &rest[1..];
        if let Some(dash) = rest.find('-') {
            let (a, b) = rest.split_at(dash);
            let b = &b[1..];
            (
                name.to_string(),
                a.parse().ok(),
                b.parse().ok(),
            )
        } else {
            (name.to_string(), rest.parse().ok(), None)
        }
    } else {
        (s.to_string(), None, None)
    }
}
