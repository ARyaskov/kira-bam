use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;

use crate::cli::FastqArgs;
use crate::collate::NameBuckets;
use crate::io::{BamReader, install_thread_pool};

const FLAG_UNMAP: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_READ1: u16 = 0x40;
const FLAG_READ2: u16 = 0x80;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

struct Out {
    w: Box<dyn Write>,
}

impl Out {
    fn new(p: &Path) -> Result<Self> {
        let f = File::create(p).with_context(|| format!("create {}", p.display()))?;
        Ok(Out {
            w: Box::new(BufWriter::with_capacity(1 << 20, f)),
        })
    }
    fn stdout() -> Self {
        Out {
            w: Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())),
        }
    }
}

pub fn run(args: FastqArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let mut r1 = args.read1.as_deref().map(Out::new).transpose()?;
    let mut r2 = args.read2.as_deref().map(Out::new).transpose()?;
    let mut singletons = args.singletons.as_deref().map(Out::new).transpose()?;
    let mut orphans = args.orphans.as_deref().map(Out::new).transpose()?;
    let mut stdout_sink = Out::stdout();

    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;

    if !args.collate {
        // Adjacent-pair scan. It needs name-grouped input, which a
        // coordinate-sorted BAM is not — every read then lands in singletons,
        // so say so rather than let an empty R1/R2 look like a successful run.
        warn_if_not_name_grouped(reader.header());

        let mut prev: Option<RecordBuf> = None;
        let mut rec = RecordBuf::default();
        while reader.read_record_buf(&mut rec)? {
            let f = u16::from(rec.flags());
            if (f & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)) != 0 {
                continue;
            }
            // The reader refills `rec` each round, so it can be moved out
            // instead of cloned.
            let cur = std::mem::take(&mut rec);
            match prev.take() {
                Some(p) => {
                    if same_qname(&p, &cur) {
                        emit_pair(
                            &p,
                            &cur,
                            &mut r1,
                            &mut r2,
                            &mut orphans,
                            &mut stdout_sink,
                            &args,
                        )?;
                    } else {
                        emit_single(&p, &mut singletons, &mut r1, &mut stdout_sink, &args)?;
                        prev = Some(cur);
                    }
                }
                None => prev = Some(cur),
            }
        }
        if let Some(p) = prev.take() {
            emit_single(&p, &mut singletons, &mut r1, &mut stdout_sink, &args)?;
        }
    } else {
        // Group by name through the collate buckets rather than sorting the
        // whole file in memory, which needed the entire BAM resident (cloned,
        // at that) before a single read could be written.
        let header = reader.header().clone();
        let mut groups = NameBuckets::new(args.tmpdir.as_deref(), &header, args.buckets.max(1))?;
        let mut rec = RecordBuf::default();
        while reader.read_record_buf(&mut rec)? {
            let f = u16::from(rec.flags());
            if (f & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)) != 0 {
                continue;
            }
            groups.push(&rec)?;
        }
        groups.for_each_group(|recs| {
            if recs.len() == 2 {
                emit_pair(
                    &recs[0],
                    &recs[1],
                    &mut r1,
                    &mut r2,
                    &mut orphans,
                    &mut stdout_sink,
                    &args,
                )
            } else {
                for rec in recs {
                    emit_single(rec, &mut singletons, &mut r1, &mut stdout_sink, &args)?;
                }
                Ok(())
            }
        })?;
    }

    Ok(())
}

/// Warn when the header does not claim name-grouped order, which the
/// adjacent-pair scan requires.
fn warn_if_not_name_grouped(header: &noodles_sam::Header) {
    use noodles_sam::header::record::value::map::header::tag as hd_tag;

    let hd = header.header();
    let value = |t| {
        hd.and_then(|h| h.other_fields().get(&t))
            .map(|v| v.to_vec())
            .unwrap_or_default()
    };
    let so = value(hd_tag::SORT_ORDER);
    let go = value(hd_tag::GROUP_ORDER);
    let grouped = so == b"queryname" || go == b"query";
    if !grouped {
        let shown = if so.is_empty() {
            "unspecified".to_string()
        } else {
            String::from_utf8_lossy(&so).into_owned()
        };
        eprintln!(
            "[kira-bam fastq] warning: @HD sort order is {shown}, not queryname. Mates are only \
             paired when they are adjacent, so this run will send most reads to singletons — \
             pass --collate, or name-sort the input first."
        );
    }
}

fn same_qname(a: &RecordBuf, b: &RecordBuf) -> bool {
    let an: Option<&[u8]> = a.name().map(AsRef::as_ref);
    let bn: Option<&[u8]> = b.name().map(AsRef::as_ref);
    an == bn
}

fn emit_pair(
    a: &RecordBuf,
    b: &RecordBuf,
    r1: &mut Option<Out>,
    r2: &mut Option<Out>,
    orphans: &mut Option<Out>,
    stdout_sink: &mut Out,
    args: &FastqArgs,
) -> Result<()> {
    let a_f = u16::from(a.flags());
    let b_f = u16::from(b.flags());
    let (first, second) = if (a_f & FLAG_READ1) != 0 {
        (a, b)
    } else if (b_f & FLAG_READ1) != 0 {
        (b, a)
    } else {
        // Unflagged pair → orphans or stdout.
        let sink = orphans.as_mut().unwrap_or(stdout_sink);
        emit_record(a, sink, args, None)?;
        emit_record(b, sink, args, None)?;
        return Ok(());
    };
    match (r1.as_mut(), r2.as_mut()) {
        (Some(out1), Some(out2)) => {
            emit_record(first, out1, args, Some(1))?;
            emit_record(second, out2, args, Some(2))?;
        }
        _ => {
            let sink = orphans.as_mut().unwrap_or(stdout_sink);
            emit_record(first, sink, args, Some(1))?;
            emit_record(second, sink, args, Some(2))?;
        }
    }
    Ok(())
}

fn emit_single(
    rec: &RecordBuf,
    singletons: &mut Option<Out>,
    r1: &mut Option<Out>,
    stdout_sink: &mut Out,
    args: &FastqArgs,
) -> Result<()> {
    let f = u16::from(rec.flags());
    let mate_num = if (f & FLAG_READ1) != 0 {
        Some(1)
    } else if (f & FLAG_READ2) != 0 {
        Some(2)
    } else {
        None
    };
    let sink = singletons.as_mut().or(r1.as_mut()).unwrap_or(stdout_sink);
    emit_record(rec, sink, args, mate_num)
}

fn emit_record(
    rec: &RecordBuf,
    out: &mut Out,
    args: &FastqArgs,
    mate_num: Option<u8>,
) -> Result<()> {
    let name: &[u8] = rec.name().map(AsRef::as_ref).unwrap_or(&[]);
    let mut suffix = String::new();
    if args.casava
        && let Some(n) = mate_num
    {
        suffix = format!("/{n}");
    }
    let flags = u16::from(rec.flags());
    if (flags & FLAG_UNMAP) == 0 && (flags & FLAG_REVERSE) != 0 {
        let seq: Vec<u8> = rec
            .sequence()
            .as_ref()
            .iter()
            .rev()
            .map(|&b| complement(b))
            .collect();
        let qual = phred33(rec, seq.len());
        write_seq(out, name, &suffix, &seq, &qual, args.fasta)?;
    } else {
        let seq: Vec<u8> = rec.sequence().as_ref().to_vec();
        let qual = phred33(rec, seq.len());
        write_seq(out, name, &suffix, &seq, &qual, args.fasta)?;
    }
    Ok(())
}

/// Phred+33 line for the record, reversed to match a reverse-complemented
/// sequence.
///
/// A record with no qualities stores 0xFF per base, which `q + 33` overflowed
/// — a panic in debug, a run of `' '` in release. FASTQ spells "unknown" as
/// `!` (Phred 0), and the line has to be exactly as long as SEQ or the record
/// is malformed, which is why the length comes from the caller.
fn phred33(rec: &RecordBuf, seq_len: usize) -> Vec<u8> {
    const MISSING: u8 = 0xFF;
    let quals = rec.quality_scores().as_ref();
    let reverse =
        (u16::from(rec.flags()) & FLAG_REVERSE) != 0 && (u16::from(rec.flags()) & FLAG_UNMAP) == 0;
    if quals.is_empty() || quals.iter().all(|&q| q == MISSING) {
        return vec![b'!'; seq_len];
    }
    let mut out: Vec<u8> = quals
        .iter()
        .map(|&q| {
            if q == MISSING {
                b'!'
            } else {
                q.saturating_add(33).min(126)
            }
        })
        .collect();
    if reverse {
        out.reverse();
    }
    out.resize(seq_len, b'!');
    out
}

fn write_seq(
    out: &mut Out,
    name: &[u8],
    suffix: &str,
    seq: &[u8],
    qual: &[u8],
    fasta: bool,
) -> Result<()> {
    if fasta {
        out.w.write_all(b">")?;
        out.w.write_all(name)?;
        out.w.write_all(suffix.as_bytes())?;
        out.w.write_all(b"\n")?;
        out.w.write_all(seq)?;
        out.w.write_all(b"\n")?;
    } else {
        out.w.write_all(b"@")?;
        out.w.write_all(name)?;
        out.w.write_all(suffix.as_bytes())?;
        out.w.write_all(b"\n")?;
        out.w.write_all(seq)?;
        out.w.write_all(b"\n+\n")?;
        out.w.write_all(qual)?;
        out.w.write_all(b"\n")?;
    }
    Ok(())
}

/// Full IUPAC complement. Passing the ambiguity codes through unchanged, as
/// this used to, silently produced a wrong reverse complement for any read
/// carrying one — BAM's 4-bit encoding represents all sixteen.
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'U' => b'A',
        b'R' => b'Y', // A/G  -> T/C
        b'Y' => b'R',
        b'S' => b'S', // C/G is its own complement
        b'W' => b'W', // A/T likewise
        b'K' => b'M', // G/T  -> C/A
        b'M' => b'K',
        b'B' => b'V', // not-A -> not-T
        b'V' => b'B',
        b'D' => b'H', // not-C -> not-G
        b'H' => b'D',
        b'N' => b'N',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'u' => b'a',
        b'r' => b'y',
        b'y' => b'r',
        b's' => b's',
        b'w' => b'w',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        b'n' => b'n',
        x => x,
    }
}
