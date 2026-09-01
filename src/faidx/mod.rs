use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

use anyhow::{Context, Result};

use crate::cli::FaidxArgs;

#[derive(Debug, Clone)]
pub(crate) struct FaiEntry {
    name: String,
    length: u64,
    offset: u64,
    line_blen: u64,
    line_len: u64,
}

pub fn run(args: FaidxArgs) -> Result<()> {
    let fai = load_or_build_fai(&args.fasta)?;

    if args.regions.is_empty() {
        return Ok(());
    }

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    let mut fasta_f = File::open(&args.fasta)?;
    let width = args.line_width.max(1);
    let mut buf: Vec<u8> = Vec::new();

    for region in &args.regions {
        let (name, start, end) = parse_region(region);
        let entry = match fai.iter().find(|e| e.name == name) {
            Some(e) => e,
            None => {
                eprintln!("[kira-bam faidx] reference {name} not found");
                continue;
            }
        };

        // Clamp to the sequence, as samtools does. Unclamped, `chr:0-100`
        // underflowed `start - 1` and a too-large end read line breaks and the
        // next record's header straight into the output.
        let start = start.unwrap_or(1).max(1);
        let end = end.unwrap_or(entry.length).min(entry.length);
        writeln!(out, ">{name}:{start}-{end}")?;
        if start > end {
            eprintln!(
                "[kira-bam faidx] {region}: empty range after clamping to {name}:1-{} \
                 (start {start} > end {end})",
                entry.length
            );
            continue;
        }

        // Wrap at `width` across the whole region rather than per file chunk:
        // the reads are bounded by the FASTA's own line length, so wrapping
        // inside the loop emitted short lines wherever those disagreed.
        let mut col = 0usize;
        let mut remaining = (end - start + 1) as usize;
        let mut pos_in_seq = start - 1;
        while remaining > 0 {
            let chunk = read_at(&mut fasta_f, entry, pos_in_seq, remaining, &mut buf)?;
            if chunk == 0 {
                break;
            }
            let mut rest: &[u8] = &buf[..chunk];
            while !rest.is_empty() {
                let take = (width - col).min(rest.len());
                out.write_all(&rest[..take])?;
                col += take;
                rest = &rest[take..];
                if col == width {
                    out.write_all(b"\n")?;
                    col = 0;
                }
            }
            pos_in_seq += chunk as u64;
            remaining -= chunk;
        }
        if col > 0 {
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Read up to one FASTA line's worth of bases starting at `pos_in_seq`
/// (0-based within the sequence) into `buf`, returning how many landed there.
fn read_at(
    f: &mut File,
    entry: &FaiEntry,
    pos_in_seq: u64,
    remaining: usize,
    buf: &mut Vec<u8>,
) -> Result<usize> {
    if entry.line_blen == 0 {
        anyhow::bail!(
            "{}: .fai reports a line length of 0; the index is malformed",
            entry.name
        );
    }
    let line_idx = pos_in_seq / entry.line_blen;
    let in_line = pos_in_seq % entry.line_blen;
    let file_offset = entry.offset + line_idx * entry.line_len + in_line;
    let chunk = (entry.line_blen - in_line).min(remaining as u64) as usize;
    buf.resize(chunk, 0);
    f.seek(SeekFrom::Start(file_offset))?;
    f.read_exact(&mut buf[..chunk])
        .with_context(|| format!("read {} at offset {file_offset}", entry.name))?;
    Ok(chunk)
}

/// Load the `.fai`, building it when absent — or when it predates the FASTA.
/// A stale index still parses, so its wrong offsets would otherwise be used
/// silently and return the wrong bases.
pub(crate) fn load_or_build_fai(fasta: &std::path::Path) -> Result<Vec<FaiEntry>> {
    let fai_path = build_fai_path(fasta);
    if fai_path.exists() {
        if !is_stale(fasta, &fai_path) {
            return read_fai(&fai_path).context("read .fai");
        }
        eprintln!(
            "[kira-bam faidx] {} is older than {} — rebuilding",
            fai_path.display(),
            fasta.display()
        );
    }
    let entries = build_fai(fasta).context("build .fai")?;
    if let Err(e) = write_fai(&fai_path, &entries) {
        // An unwritable directory is no reason to fail the query itself.
        eprintln!(
            "[kira-bam faidx] could not write {}: {e:#}",
            fai_path.display()
        );
    }
    Ok(entries)
}

/// Was `derived` last modified before `source`? Unknown timestamps count as
/// fresh — a missing mtime is not evidence of staleness.
pub(crate) fn is_stale(source: &std::path::Path, derived: &std::path::Path) -> bool {
    let modified = |p: &std::path::Path| std::fs::metadata(p).and_then(|m| m.modified()).ok();
    match (modified(source), modified(derived)) {
        (Some(src), Some(der)) => der < src,
        _ => false,
    }
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
    // Split on the *last* colon and only if the tail parses as a range: contig
    // names may contain colons themselves (HLA alleles, for one).
    let Some((name, rest)) = s.rsplit_once(':') else {
        return (s.to_string(), None, None);
    };
    let rest = rest.replace(',', "");
    let parsed = match rest.split_once('-') {
        Some((a, b)) => {
            let start = if a.is_empty() { None } else { a.parse().ok() };
            let end = if b.is_empty() { None } else { b.parse().ok() };
            let start_ok = a.is_empty() || start.is_some();
            let end_ok = b.is_empty() || end.is_some();
            if start_ok && end_ok && !(a.is_empty() && b.is_empty()) {
                Some((start, end))
            } else {
                None
            }
        }
        None => rest.parse().ok().map(|s| (Some(s), None)),
    };
    match parsed {
        Some((start, end)) => (name.to_string(), start, end),
        None => (s.to_string(), None, None),
    }
}

/// Random-access FASTA backed by its `.fai`, caching one contig at a time.
///
/// Coordinate-sorted input touches references in header order, so holding only
/// the current one bounds memory by the largest contig (~250 MB for human
/// chr1) instead of the whole genome (~3.1 GB, which loading every record into
/// a map cost — twice over, since uppercasing copied it again).
pub struct RefCache {
    file: File,
    entries: Vec<FaiEntry>,
    /// Index into `entries` of the contig currently in `seq`.
    cached: Option<usize>,
    /// Uppercased bases of the cached contig, newlines stripped.
    seq: Vec<u8>,
    /// Raw file bytes, reused across loads.
    raw: Vec<u8>,
}

impl RefCache {
    pub fn open(fasta: &std::path::Path) -> Result<Self> {
        let entries = load_or_build_fai(fasta)?;
        let file =
            File::open(fasta).with_context(|| format!("open reference {}", fasta.display()))?;
        Ok(Self {
            file,
            entries,
            cached: None,
            seq: Vec::new(),
            raw: Vec::new(),
        })
    }

    /// The full sequence of `name`, loading it if another contig is cached.
    /// `None` when the FASTA has no such record.
    pub fn sequence(&mut self, name: &str) -> Result<Option<&[u8]>> {
        let Some(idx) = self.entries.iter().position(|e| e.name == name) else {
            return Ok(None);
        };
        if self.cached != Some(idx) {
            self.load(idx)?;
            self.cached = Some(idx);
        }
        Ok(Some(&self.seq))
    }

    fn load(&mut self, idx: usize) -> Result<()> {
        let FaiEntry {
            length,
            offset,
            line_blen,
            line_len,
            ..
        } = self.entries[idx];
        let name = self.entries[idx].name.clone();

        self.seq.clear();
        if length == 0 {
            return Ok(());
        }
        if line_blen == 0 {
            anyhow::bail!("{name}: .fai reports a line length of 0; the index is malformed");
        }
        // One contiguous read of the record's byte span beats a seek per line;
        // the line terminators come out in the copy below.
        let full_lines = length / line_blen;
        let remainder = length % line_blen;
        let mut span = full_lines * line_len + remainder;
        if remainder == 0 {
            // The final line's terminator is not part of the sequence.
            span -= line_len - line_blen;
        }

        let mut raw = std::mem::take(&mut self.raw);
        raw.resize(span as usize, 0);
        let read = self
            .file
            .seek(SeekFrom::Start(offset))
            .and_then(|_| self.file.read_exact(&mut raw));
        if let Err(e) = read {
            self.raw = raw;
            return Err(anyhow::Error::new(e))
                .with_context(|| format!("read {name} from the reference"));
        }

        self.seq.reserve(length as usize);
        self.seq.extend(
            raw.iter()
                .filter(|b| !matches!(b, b'\n' | b'\r'))
                .map(|b| b.to_ascii_uppercase()),
        );
        self.raw = raw;
        Ok(())
    }
}
