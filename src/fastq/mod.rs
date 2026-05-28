use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;

use crate::cli::FastqArgs;
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

    let mut reader = BamReader::open(&args.input).context("open input")?;

    if !args.collate {
        let mut prev: Option<RecordBuf> = None;
        let mut rec = RecordBuf::default();
        while reader.read_record_buf(&mut rec)? {
            let f = u16::from(rec.flags());
            if (f & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)) != 0 {
                continue;
            }
            match prev.take() {
                Some(p) => {
                    if same_qname(&p, &rec) {
                        emit_pair(
                            &p,
                            &rec,
                            &mut r1,
                            &mut r2,
                            &mut orphans,
                            &mut stdout_sink,
                            &args,
                        )?;
                    } else {
                        emit_single(&p, &mut singletons, &mut r1, &mut stdout_sink, &args)?;
                        prev = Some(rec.clone());
                    }
                }
                None => prev = Some(rec.clone()),
            }
        }
        if let Some(p) = prev.take() {
            emit_single(&p, &mut singletons, &mut r1, &mut stdout_sink, &args)?;
        }
    } else {
        let mut buf: Vec<RecordBuf> = Vec::new();
        let mut rec = RecordBuf::default();
        while reader.read_record_buf(&mut rec)? {
            let f = u16::from(rec.flags());
            if (f & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)) != 0 {
                continue;
            }
            buf.push(rec.clone());
        }
        buf.sort_by(|a, b| {
            let an: &[u8] = a.name().map(AsRef::as_ref).unwrap_or(&[]);
            let bn: &[u8] = b.name().map(AsRef::as_ref).unwrap_or(&[]);
            an.cmp(bn)
        });
        let mut i = 0;
        while i < buf.len() {
            let j = (i + 1..buf.len())
                .find(|&k| !same_qname(&buf[i], &buf[k]))
                .unwrap_or(buf.len());
            if j - i == 2 {
                let (a, b) = (buf[i].clone(), buf[i + 1].clone());
                emit_pair(
                    &a,
                    &b,
                    &mut r1,
                    &mut r2,
                    &mut orphans,
                    &mut stdout_sink,
                    &args,
                )?;
            } else {
                for rec in buf.iter().take(j).skip(i) {
                    let rec = rec.clone();
                    emit_single(&rec, &mut singletons, &mut r1, &mut stdout_sink, &args)?;
                }
            }
            i = j;
        }
    }

    Ok(())
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
        let qual: Vec<u8> = rec
            .quality_scores()
            .as_ref()
            .iter()
            .rev()
            .map(|&q| q + 33)
            .collect();
        write_seq(out, name, &suffix, &seq, &qual, args.fasta)?;
    } else {
        let seq: Vec<u8> = rec.sequence().as_ref().to_vec();
        let qual: Vec<u8> = rec
            .quality_scores()
            .as_ref()
            .iter()
            .map(|&q| q + 33)
            .collect();
        write_seq(out, name, &suffix, &seq, &qual, args.fasta)?;
    }
    Ok(())
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

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'N' | b'n' => b'N',
        x => x,
    }
}
