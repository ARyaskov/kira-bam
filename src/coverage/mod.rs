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

#[derive(Default)]
struct RefStats {
    name: String,
    length: usize,
    numreads: u64,
    /// Bitmap (one u8 per base, 1 = covered). Compact at 1 byte/base; for chr1 ≈ 250 MB.
    /// Fine in practice, awful for whole-genome single-ref.
    /// TODO: use BitSet for memory.
    covered: Vec<u8>,
    sum_depth: u64,
    sum_baseq: u64,
    sum_baseq_n: u64,
    sum_mapq: u64,
    sum_mapq_n: u64,
}

pub fn run(args: CoverageArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let mut reader = BamReader::open(&args.input).context("open input")?;
    let header = reader.header().clone();

    let mut stats: HashMap<usize, RefStats, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();

    // Pre-populate so refs with zero reads still get a row.
    for (i, (name, sq)) in header.reference_sequences().iter().enumerate() {
        let len = usize::from(sq.length());
        stats.insert(
            i,
            RefStats {
                name: name.to_string(),
                length: len,
                covered: vec![0u8; len],
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
        if let Some(mapq) = rec.mapping_quality() {
            if u8::from(mapq) < args.min_mapq {
                continue;
            }
        }
        let tid = match rec.reference_sequence_id() {
            Some(t) => t,
            None => continue,
        };
        let start = match rec.alignment_start() {
            Some(p) => usize::from(p),
            None => continue,
        };
        let end = match rec.alignment_end() {
            Some(p) => usize::from(p),
            None => continue,
        };
        let st = stats.get_mut(&tid).unwrap();
        st.numreads += 1;
        let mapq = rec.mapping_quality().map(u8::from).unwrap_or(0) as u64;
        st.sum_mapq += mapq;
        st.sum_mapq_n += 1;

        // Per-base mark + baseq sum (approximate: just count bases by quality threshold).
        let qs = rec.quality_scores().as_ref();
        let mut q_iter = qs.iter().copied();
        for pos in start..=end.min(st.length) {
            if pos == 0 || pos > st.length {
                continue;
            }
            st.covered[pos - 1] = 1;
            st.sum_depth += 1;
            if let Some(q) = q_iter.next() {
                if q >= args.min_baseq {
                    st.sum_baseq += q as u64;
                    st.sum_baseq_n += 1;
                }
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
        if !args.regions.is_empty() && !args.regions.iter().any(|r| r == &st.name) {
            continue;
        }
        let covbases: u64 = st.covered.iter().map(|&b| b as u64).sum();
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
