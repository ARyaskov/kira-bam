use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};

use anyhow::{Context, Result};
use rustc_hash::FxHasher;

use crate::cli::FqidxArgs;
use crate::fastq_input::FastqInput;
use crate::io::BufferedOutput;

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

    let by_name: HashMap<String, &FqiEntry, std::hash::BuildHasherDefault<FxHasher>> =
        entries.iter().map(|e| (e.name.clone(), e)).collect();

    let mut out = BufferedOutput::create(args.output.as_deref())?;
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
    out.finish()
}

fn build_fqi_path(fq: &std::path::Path) -> std::path::PathBuf {
    let mut s = fq.as_os_str().to_os_string();
    s.push(".fqi");
    std::path::PathBuf::from(s)
}

fn build_fqi(path: &std::path::Path) -> Result<Vec<FqiEntry>> {
    let mut reader = FastqInput::open(path)?;
    let mut entries = Vec::new();
    while let Some(record) = reader.next()? {
        entries.push(FqiEntry {
            name: record.name.to_owned(),
            length: record.sequence.len() as u64,
            seq_offset: record.sequence_offset,
            qual_offset: record.quality_offset,
        });
    }
    Ok(entries)
}

fn write_fqi(path: &std::path::Path, entries: &[FqiEntry]) -> Result<()> {
    let mut f = BufferedOutput::create(Some(path))?;
    for e in entries {
        writeln!(
            f,
            "{}\t{}\t{}\t{}",
            e.name, e.length, e.seq_offset, e.qual_offset
        )?;
    }
    f.finish()
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
