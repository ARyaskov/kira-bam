use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record::data::field::Tag;

use crate::cli::FixmateArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAP: u16 = 0x4;
const FLAG_MUNMAP: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

pub fn run(args: FixmateArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let mut reader = BamReader::open(&args.input).context("open input")?;
    let mut header = reader.header().clone();
    append_pg(&mut header, &PgInfo::new("fixmate", !args.no_pg))?;

    let fmt = if args.uncompressed {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    };
    let mut writer = BamWriter::create(Some(&args.output), header, fmt)?;
    writer.write_header()?;

    let mut current_qname: Vec<u8> = Vec::new();
    let mut group: Vec<RecordBuf> = Vec::new();
    let mut rec = RecordBuf::default();

    while reader.read_record_buf(&mut rec)? {
        let q: Vec<u8> = rec.name().map(|n| n.to_vec()).unwrap_or_default();
        if q != current_qname && !group.is_empty() {
            flush_group(&mut group, &mut writer, &args)?;
            current_qname.clear();
        }
        if current_qname.is_empty() {
            current_qname = q;
        }
        group.push(rec.clone());
    }
    if !group.is_empty() {
        flush_group(&mut group, &mut writer, &args)?;
    }
    writer.finish()?;
    Ok(())
}

fn flush_group(
    group: &mut Vec<RecordBuf>,
    writer: &mut BamWriter,
    args: &FixmateArgs,
) -> Result<()> {
    let primary: Vec<usize> = (0..group.len())
        .filter(|&i| {
            let f = u16::from(group[i].flags());
            (f & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)) == 0
        })
        .collect();

    if primary.len() == 2 {
        let (a_idx, b_idx) = (primary[0], primary[1]);
        let (left, right) = group.split_at_mut(b_idx);
        let a = &mut left[a_idx];
        let b = &mut right[0];
        fix_mate_pair(a, b, args);
    } else if primary.len() == 1 {
        // Singleton: mark as unpaired-mate.
        let r = &mut group[primary[0]];
        let mut f = u16::from(r.flags());
        f |= FLAG_MUNMAP;
        f &= !FLAG_PROPER_PAIR;
        *r.flags_mut() = Flags::from(f);
    }

    for r in group.drain(..) {
        let f = u16::from(r.flags());
        if args.remove_unpaired {
            if (f & FLAG_UNMAP) != 0 {
                continue;
            }
            if (f & FLAG_PAIRED) != 0 && (f & FLAG_MUNMAP) != 0 {
                continue;
            }
        }
        let _ = writer.write_record_buf(&r);
    }
    Ok(())
}

fn fix_mate_pair(a: &mut RecordBuf, b: &mut RecordBuf, args: &FixmateArgs) {
    let a_flags = u16::from(a.flags());
    let b_flags = u16::from(b.flags());
    let a_unmap = (a_flags & FLAG_UNMAP) != 0;
    let b_unmap = (b_flags & FLAG_UNMAP) != 0;

    // Cross-set mate reference/position.
    *a.mate_reference_sequence_id_mut() = b.reference_sequence_id();
    *a.mate_alignment_start_mut() = b.alignment_start();
    *b.mate_reference_sequence_id_mut() = a.reference_sequence_id();
    *b.mate_alignment_start_mut() = a.alignment_start();

    // Cross-set mate strand and mate-unmap bits.
    let mut na = a_flags;
    let mut nb = b_flags;
    if (b_flags & FLAG_REVERSE) != 0 {
        na |= FLAG_MATE_REVERSE;
    } else {
        na &= !FLAG_MATE_REVERSE;
    }
    if (a_flags & FLAG_REVERSE) != 0 {
        nb |= FLAG_MATE_REVERSE;
    } else {
        nb &= !FLAG_MATE_REVERSE;
    }
    if b_unmap {
        na |= FLAG_MUNMAP;
    } else {
        na &= !FLAG_MUNMAP;
    }
    if a_unmap {
        nb |= FLAG_MUNMAP;
    } else {
        nb &= !FLAG_MUNMAP;
    }

    // Proper-pair: same ref, both mapped, opposite strands, reasonable distance.
    let proper = !a_unmap
        && !b_unmap
        && a.reference_sequence_id() == b.reference_sequence_id()
        && a.reference_sequence_id().is_some()
        && ((a_flags ^ b_flags) & FLAG_REVERSE) != 0;
    if proper {
        na |= FLAG_PROPER_PAIR;
        nb |= FLAG_PROPER_PAIR;
    } else {
        na &= !FLAG_PROPER_PAIR;
        nb &= !FLAG_PROPER_PAIR;
    }

    *a.flags_mut() = Flags::from(na);
    *b.flags_mut() = Flags::from(nb);

    // TLEN: signed insert size.
    if a.reference_sequence_id() == b.reference_sequence_id() && !a_unmap && !b_unmap {
        let a_start = a.alignment_start().map(|p| usize::from(p) as i64).unwrap_or(0);
        let b_start = b.alignment_start().map(|p| usize::from(p) as i64).unwrap_or(0);
        let a_end = a.alignment_end().map(|p| usize::from(p) as i64).unwrap_or(a_start);
        let b_end = b.alignment_end().map(|p| usize::from(p) as i64).unwrap_or(b_start);
        let (left_start, right_end) = if a_start <= b_start {
            (a_start, b_end.max(a_end))
        } else {
            (b_start, a_end.max(b_end))
        };
        let isize_val = right_end - left_start + 1;
        let a_tlen: i32 = if a_start <= b_start {
            isize_val as i32
        } else {
            -(isize_val as i32)
        };
        *a.template_length_mut() = a_tlen;
        *b.template_length_mut() = -a_tlen;
    } else {
        *a.template_length_mut() = 0;
        *b.template_length_mut() = 0;
    }

    if args.add_mc_ms {
        add_mc_ms_tags(a, b);
    }
}

fn add_mc_ms_tags(a: &mut RecordBuf, b: &mut RecordBuf) {
    // MC = mate CIGAR string
    let a_cigar = format_cigar(a);
    let b_cigar = format_cigar(b);
    if !b_cigar.is_empty() {
        a.data_mut()
            .insert(Tag::MATE_CIGAR, Value::String(b_cigar.into()));
    }
    if !a_cigar.is_empty() {
        b.data_mut()
            .insert(Tag::MATE_CIGAR, Value::String(a_cigar.into()));
    }
    // MS = sum of mate base quality scores (capped at >=15 to match samtools convention)
    let a_ms = quality_sum_filtered(a);
    let b_ms = quality_sum_filtered(b);
    a.data_mut()
        .insert(Tag::new(b'm', b's'), Value::Int32(b_ms as i32));
    b.data_mut()
        .insert(Tag::new(b'm', b's'), Value::Int32(a_ms as i32));
}

fn format_cigar(rec: &RecordBuf) -> String {
    use noodles_sam::alignment::record::cigar::op::Kind;
    let mut s = String::new();
    for op in rec.cigar().as_ref().iter() {
        let c = match op.kind() {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        };
        use std::fmt::Write;
        let _ = write!(s, "{}{}", op.len(), c);
    }
    s
}

fn quality_sum_filtered(rec: &RecordBuf) -> u32 {
    rec.quality_scores()
        .as_ref()
        .iter()
        .filter(|&&q| q >= 15)
        .map(|&q| q as u32)
        .sum()
}
