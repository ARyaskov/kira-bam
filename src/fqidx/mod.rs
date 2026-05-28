use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

use anyhow::{Context, Result};
use rustc_hash::FxHasher;

use crate::cli::FqidxArgs;

/// `.fqi` per-record: name \t length \t seq_offset \t qual_offset
#[derive(Debug, Clone)]
struct FqiEntry {
    name: String,
    length: u64,
    seq_offset: u64,
    qual_offset: u64,
}

pub fn run(args: FqidxArgs) -> Result<()> {
    let fqi_path = build_fqi_path(&args.fastq);
    let entries: Vec<FqiEntry> = if !fqi_path.exists() {
        let e = build_fqi(&args.fastq).context("build .fqi")?;
        write_fqi(&fqi_path, &e).context("write .fqi")?;
        e
    } else {
        read_fqi(&fqi_path).context("read .fqi")?
    };

    if args.names.is_empty() {
        return Ok(());
    }

    let by_name: HashMap<String, &FqiEntry, std::hash::BuildHasherDefault<FxHasher>> = entries
        .iter()
        .map(|e| (e.name.clone(), e))
        .collect();

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    let mut f = File::open(&args.fastq)?;

    for name in &args.names {
        let entry = match by_name.get(name.as_str()) {
            Some(e) => *e,
            None => {
                eprintln!("[kira-bam fqidx] not found: {name}");
                continue;
            }
        };
        // Seek + read sequence
        f.seek(SeekFrom::Start(entry.seq_offset))?;
        let mut seq = vec![0u8; entry.length as usize];
        f.read_exact(&mut seq)?;
        f.seek(SeekFrom::Start(entry.qual_offset))?;
        let mut qual = vec![0u8; entry.length as usize];
        f.read_exact(&mut qual)?;

        writeln!(out, "@{}", entry.name)?;
        out.write_all(&seq)?;
        out.write_all(b"\n+\n")?;
        out.write_all(&qual)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn build_fqi_path(fq: &std::path::Path) -> std::path::PathBuf {
    let mut s = fq.as_os_str().to_os_string();
    s.push(".fqi");
    std::path::PathBuf::from(s)
}

fn build_fqi(path: &std::path::Path) -> Result<Vec<FqiEntry>> {
    let f = File::open(path).context("open FASTQ")?;
    let mut r = BufReader::with_capacity(1 << 20, f);

    let mut out = Vec::new();
    let mut buf = String::new();
    let mut pos: u64 = 0;
    loop {
        buf.clear();
        let n = r.read_line(&mut buf)?;
        if n == 0 {
            break;
        }
        if !buf.starts_with('@') {
            // Resilient skip — sometimes blank lines or stray content.
            pos += n as u64;
            continue;
        }
        let name = buf
            .trim_end_matches(['\n', '\r'])
            .strip_prefix('@')
            .unwrap_or("")
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        pos += n as u64;

        let seq_offset = pos;
        buf.clear();
        let n2 = r.read_line(&mut buf)?;
        let seq_len = buf.trim_end_matches(['\n', '\r']).len() as u64;
        pos += n2 as u64;

        // `+` line
        buf.clear();
        let n3 = r.read_line(&mut buf)?;
        pos += n3 as u64;

        let qual_offset = pos;
        buf.clear();
        let n4 = r.read_line(&mut buf)?;
        pos += n4 as u64;

        out.push(FqiEntry {
            name,
            length: seq_len,
            seq_offset,
            qual_offset,
        });
    }
    Ok(out)
}

fn write_fqi(path: &std::path::Path, entries: &[FqiEntry]) -> Result<()> {
    let mut f = BufWriter::new(File::create(path)?);
    for e in entries {
        writeln!(
            f,
            "{}\t{}\t{}\t{}",
            e.name, e.length, e.seq_offset, e.qual_offset
        )?;
    }
    Ok(())
}

fn read_fqi(path: &std::path::Path) -> Result<Vec<FqiEntry>> {
    let f = File::open(path)?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    for line in r.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }
        out.push(FqiEntry {
            name: parts[0].to_string(),
            length: parts[1].parse().unwrap_or(0),
            seq_offset: parts[2].parse().unwrap_or(0),
            qual_offset: parts[3].parse().unwrap_or(0),
        });
    }
    Ok(out)
}
