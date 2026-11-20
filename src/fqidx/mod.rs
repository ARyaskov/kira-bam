use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};

use anyhow::{Context, Result, ensure};
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
    let mut f =
        File::open(&args.fastq).with_context(|| format!("open FASTQ {}", args.fastq.display()))?;
    let entries = if !fqi_path.try_exists()? || crate::faidx::is_stale(&args.fastq, &fqi_path) {
        let e = build_fqi(&args.fastq).context("build .fqi")?;
        write_fqi(&fqi_path, &e).context("write .fqi")?;
        e
    } else {
        read_fqi(&fqi_path, f.metadata()?.len())?
    };

    if args.names.is_empty() {
        return Ok(());
    }

    let by_name: HashMap<String, &FqiEntry, std::hash::BuildHasherDefault<FxHasher>> =
        entries.iter().map(|e| (e.name.clone(), e)).collect();

    let mut out = BufferedOutput::create(args.output.as_deref())?;

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
        let length = usize::try_from(entry.length).context("FASTQ read length exceeds usize")?;
        let mut seq = vec![0u8; length];
        f.read_exact(&mut seq)?;
        f.seek(SeekFrom::Start(entry.qual_offset))?;
        let mut qual = vec![0u8; length];
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

fn read_fqi(path: &std::path::Path, fastq_len: u64) -> Result<Vec<FqiEntry>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    for (index, line) in r.lines().enumerate() {
        let location = || format!("{}: line {}", path.display(), index + 1);
        let line = line.with_context(location)?;
        out.push(parse_entry(&line, fastq_len).with_context(location)?);
    }
    Ok(out)
}

fn parse_entry(line: &str, fastq_len: u64) -> Result<FqiEntry> {
    let fields: Vec<_> = line.split('\t').collect();
    let [name, length, seq_offset, qual_offset] = fields.as_slice() else {
        anyhow::bail!("expected four .fqi fields");
    };
    ensure!(!name.is_empty(), "empty read name");
    let length: u64 = length.parse().context("invalid read length")?;
    let seq_offset: u64 = seq_offset.parse().context("invalid sequence offset")?;
    let qual_offset: u64 = qual_offset.parse().context("invalid quality offset")?;
    ensure!(
        seq_offset > 0
            && qual_offset
                .checked_sub(seq_offset)
                .is_some_and(|span| span > length)
            && qual_offset
                .checked_add(length)
                .is_some_and(|end| end <= fastq_len),
        "read offsets or length are inconsistent with the FASTQ file"
    );
    Ok(FqiEntry {
        name: (*name).to_owned(),
        length,
        seq_offset,
        qual_offset,
    })
}
