use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
use rustc_hash::FxHasher;
use std::collections::HashMap;

use crate::cli::CoverageArgs;
use crate::io::{BamReader, install_thread_pool};

const ALL_COLUMNS: [&str; 9] = [
    "rname",
    "startpos",
    "endpos",
    "numreads",
    "covbases",
    "coverage",
    "meandepth",
    "meanbaseq",
    "meanmapq",
];

/// One bit per reference base, packed 64 to a word.
///
/// A byte per base cost ~3.1 GB across a human reference — and the whole
/// reference was allocated up front, `-r` or not. Bits cut that eightfold, and
/// only the contigs actually reported get allocated at all.
#[derive(Default)]
struct CoveredBits {
    words: Vec<u64>,
    /// Bases set so far, counted as they are set rather than by a final scan
    /// over every base of the contig.
    popcount: u64,
}

impl CoveredBits {
    fn with_len(len: usize) -> Self {
        Self {
            words: vec![0u64; len.div_ceil(64)],
            popcount: 0,
        }
    }

    /// Mark the 0-based position `i` covered.
    fn set(&mut self, i: usize) {
        let (word, bit) = (i / 64, i % 64);
        let mask = 1u64 << bit;
        if let Some(w) = self.words.get_mut(word)
            && *w & mask == 0
        {
            *w |= mask;
            self.popcount += 1;
        }
    }
}

#[derive(Default)]
struct RefStats {
    name: String,
    length: usize,
    numreads: u64,
    covered: CoveredBits,
    sum_depth: u64,
    sum_baseq: u64,
    sum_baseq_n: u64,
    sum_mapq: u64,
    sum_mapq_n: u64,
}

pub fn run(args: CoverageArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let header = reader.header().clone();

    let mut stats: HashMap<usize, RefStats, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();

    // Pre-populate so refs with zero reads still get a row — but only the ones
    // `-r` asks for. Allocating a bitmap for every contig regardless meant a
    // single-region run still paid for the whole genome.
    for (i, (name, sq)) in header.reference_sequences().iter().enumerate() {
        let name = name.to_string();
        if !args.regions.is_empty() && !args.regions.contains(&name) {
            continue;
        }
        let len = usize::from(sq.length());
        stats.insert(
            i,
            RefStats {
                name,
                length: len,
                covered: CoveredBits::with_len(len),
                ..Default::default()
            },
        );
    }

    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        let flags = u16::from(rec.flags());
        if (flags & 0x4) != 0 {
            continue;
        }
        if let Some(mapq) = rec.mapping_quality()
            && u8::from(mapq) < args.min_mapq
        {
            continue;
        }
        let tid = match rec.reference_sequence_id() {
            Some(t) => t,
            None => continue,
        };
        let start = match rec.alignment_start() {
            Some(p) => usize::from(p),
            None => continue,
        };
        // Absent when `-r` excluded this contig.
        let Some(st) = stats.get_mut(&tid) else {
            continue;
        };
        st.numreads += 1;
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0) as u64;
        st.sum_mapq += mapq;
        st.sum_mapq_n += 1;

        // CIGAR-aware per-base accounting. Only M/=/X positions are covered, and
        // each read base's quality is matched to the reference position it
        // actually aligns to. The previous code iterated the [start,end] span and
        // pulled qualities positionally, which (a) counted deletions and N-skips
        // (introns) as covered and (b) mis-assigned base qualities across any
        // indel. -Q now excludes low-quality bases from depth too, matching
        // samtools coverage semantics.
        let quals = rec.quality_scores().as_ref();
        let mut ref_pos = start;
        let mut read_idx = 0usize;
        for op in rec.cigar().as_ref().iter() {
            let kind = op.kind();
            let len = op.len();
            let cons_ref = kind.consumes_reference();
            let cons_read = kind.consumes_read();
            if cons_ref && cons_read {
                for _ in 0..len {
                    if ref_pos >= 1 && ref_pos <= st.length {
                        let q = quals.get(read_idx).copied().unwrap_or(0);
                        if q >= args.min_baseq {
                            st.covered.set(ref_pos - 1);
                            st.sum_depth += 1;
                            st.sum_baseq += q as u64;
                            st.sum_baseq_n += 1;
                        }
                    }
                    ref_pos += 1;
                    read_idx += 1;
                }
            } else if cons_ref {
                ref_pos += len;
            } else if cons_read {
                read_idx += len;
            }
        }
    }

    let columns = parse_columns(args.columns.as_deref())?;
    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    if !args.no_header {
        let header_line = columns
            .iter()
            .map(|c| format!("#{c}"))
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(out, "{header_line}")?;
    }

    let mut rows: Vec<&RefStats> = stats.values().collect();
    rows.sort_by_key(|r| header_index(&header, &r.name));

    for st in rows {
        // Maintained as bases are marked, rather than an O(contig length) scan
        // per row at the end.
        let covbases: u64 = st.covered.popcount;
        let coverage_pct = if st.length > 0 {
            100.0 * covbases as f64 / st.length as f64
        } else {
            0.0
        };
        let meandepth = if st.length > 0 {
            st.sum_depth as f64 / st.length as f64
        } else {
            0.0
        };
        let meanbaseq = if st.sum_baseq_n > 0 {
            st.sum_baseq as f64 / st.sum_baseq_n as f64
        } else {
            0.0
        };
        let meanmapq = if st.sum_mapq_n > 0 {
            st.sum_mapq as f64 / st.sum_mapq_n as f64
        } else {
            0.0
        };

        let cells: Vec<String> = columns
            .iter()
            .map(|c| match *c {
                "rname" => st.name.clone(),
                "startpos" => "1".to_string(),
                "endpos" => st.length.to_string(),
                "numreads" => st.numreads.to_string(),
                "covbases" => covbases.to_string(),
                "coverage" => format!("{coverage_pct:.4}"),
                "meandepth" => format!("{meandepth:.4}"),
                "meanbaseq" => format!("{meanbaseq:.4}"),
                "meanmapq" => format!("{meanmapq:.4}"),
                _ => String::new(),
            })
            .collect();
        writeln!(out, "{}", cells.join("\t"))?;
    }
    out.flush()?;
    Ok(())
}

fn parse_columns(spec: Option<&str>) -> Result<Vec<&'static str>> {
    match spec {
        None => Ok(ALL_COLUMNS.to_vec()),
        Some(s) => {
            let mut out: Vec<&'static str> = Vec::new();
            for c in s.split(',') {
                let c = c.trim();
                if let Some(&matched) = ALL_COLUMNS.iter().find(|x| **x == c) {
                    out.push(matched);
                } else {
                    anyhow::bail!("unknown column {c:?}; valid: {}", ALL_COLUMNS.join(","));
                }
            }
            Ok(out)
        }
    }
}

fn header_index(h: &Header, name: &str) -> usize {
    h.reference_sequences()
        .iter()
        .position(|(n, _)| n == name)
        .unwrap_or(usize::MAX)
}
