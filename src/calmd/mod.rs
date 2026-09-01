use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;

use crate::cli::CalmdArgs;
use crate::faidx::RefCache;
use crate::io::{
    BamReader, BamWriter, PgInfo, append_pg, install_thread_pool, report_skipped,
    write_record_tolerant,
};
use crate::types::OutputFormat;

pub fn run(args: CalmdArgs) -> Result<()> {
    install_thread_pool(args.threads);

    // One contig at a time via the `.fai`, not the whole genome in a map.
    let mut refs = RefCache::open(&args.reference)?;

    let mut reader =
        BamReader::open_with_reference(&args.input, Some(&args.reference)).context("open input")?;
    let mut header = reader.header().clone();
    append_pg(&mut header, &PgInfo::new("calmd", !args.no_pg))?;

    let fmt = if args.uncompressed {
        OutputFormat::UncompressedBam
    } else if args.bam {
        OutputFormat::Bam
    } else {
        OutputFormat::Sam
    };
    let mut writer = BamWriter::create(args.output.as_deref(), header.clone(), fmt)?;
    writer.write_header()?;

    // Contig names resolved once, not rebuilt as a String per record.
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| n.to_string())
        .collect();

    let mut rec = RecordBuf::default();
    let mut skipped = 0u64;
    let mut md = String::new();
    while reader.read_record_buf(&mut rec)? {
        if let Some(tid) = rec.reference_sequence_id() {
            let ref_name = ref_names.get(tid).map(String::as_str).unwrap_or("");
            if let Some(seq) = refs.sequence(ref_name)?
                && let Some(start) = rec.alignment_start()
            {
                let start = usize::from(start);
                if let Some(nm) = compute_md_nm(&rec, seq, start, &mut md) {
                    // samtools calmd writes exactly two tags: MD (mismatch string)
                    // and NM (edit distance). It must NOT touch NH — that tag is the
                    // number of reported alignments for a multi-mapper and is owned by
                    // the aligner. The previous code wrote NH = NM, silently corrupting
                    // multi-mapper bookkeeping for any downstream RNA-seq tool.
                    rec.data_mut()
                        .insert(Tag::MISMATCHED_POSITIONS, Value::String(md.clone().into()));
                    rec.data_mut()
                        .insert(Tag::EDIT_DISTANCE, Value::Int32(nm as i32));
                }
            }
        }
        write_record_tolerant(&mut writer, &rec, "calmd", &mut skipped)?;
    }
    writer.finish()?;
    report_skipped("calmd", skipped);
    Ok(())
}

/// Build the MD string into `md` (cleared first) and return the edit distance.
/// `None` leaves `md` unspecified and means the record could not be scored.
fn compute_md_nm(
    rec: &RecordBuf,
    refseq: &[u8],
    ref_start_1: usize,
    md: &mut String,
) -> Option<u32> {
    use std::fmt::Write;
    let seq = rec.sequence().as_ref();
    if seq.is_empty() {
        return None;
    }
    md.clear();
    let mut nm: u32 = 0;
    let mut read_idx: usize = 0;
    let mut ref_idx: usize = ref_start_1.saturating_sub(1);
    let mut run_match: u32 = 0;
    for op in rec.cigar().as_ref().iter() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    if ref_idx >= refseq.len() || read_idx >= seq.len() {
                        return None;
                    }
                    let rb = refseq[ref_idx].to_ascii_uppercase();
                    let qb = seq[read_idx].to_ascii_uppercase();
                    if rb == qb {
                        run_match += 1;
                    } else {
                        write!(md, "{run_match}").ok()?;
                        run_match = 0;
                        md.push(rb as char);
                        nm += 1;
                    }
                    ref_idx += 1;
                    read_idx += 1;
                }
            }
            Kind::Insertion => {
                read_idx += len;
                nm += len as u32;
            }
            Kind::Deletion | Kind::Skip => {
                if op.kind() == Kind::Deletion {
                    write!(md, "{run_match}").ok()?;
                    run_match = 0;
                    md.push('^');
                    for _ in 0..len {
                        if ref_idx >= refseq.len() {
                            return None;
                        }
                        md.push(refseq[ref_idx].to_ascii_uppercase() as char);
                        ref_idx += 1;
                    }
                    nm += len as u32;
                } else {
                    ref_idx += len;
                }
            }
            Kind::SoftClip => {
                read_idx += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    write!(md, "{run_match}").ok()?;
    Some(nm)
}
