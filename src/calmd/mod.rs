use std::collections::HashMap;
use std::fs::File;

use anyhow::{Context, Result};
use noodles_fasta as fasta;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use rustc_hash::FxHasher;

use crate::cli::CalmdArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

pub fn run(args: CalmdArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let refs = load_reference(&args.reference)?;

    let mut reader = BamReader::open(&args.input).context("open input")?;
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

    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        if let Some(tid) = rec.reference_sequence_id() {
            let ref_name: String = header
                .reference_sequences()
                .get_index(tid)
                .map(|(n, _)| n.to_string())
                .unwrap_or_default();
            if let Some(seq) = refs.get(ref_name.as_bytes())
                && let Some(start) = rec.alignment_start()
            {
                let start = usize::from(start);
                if let Some((md, nm)) = compute_md_nm(&rec, seq, start) {
                    rec.data_mut()
                        .insert(Tag::MISMATCHED_POSITIONS, Value::String(md.into()));
                    rec.data_mut()
                        .insert(Tag::ALIGNMENT_HIT_COUNT, Value::Int32(nm as i32));
                    rec.data_mut()
                        .insert(Tag::EDIT_DISTANCE, Value::Int32(nm as i32));
                }
            }
        }
        let _ = writer.write_record_buf(&rec);
    }
    writer.finish()?;
    Ok(())
}

fn load_reference(
    path: &std::path::Path,
) -> Result<HashMap<Vec<u8>, Vec<u8>, std::hash::BuildHasherDefault<FxHasher>>> {
    let file = File::open(path).with_context(|| format!("open reference {}", path.display()))?;
    let mut reader = fasta::io::Reader::new(std::io::BufReader::with_capacity(1 << 20, file));
    let mut out: HashMap<Vec<u8>, Vec<u8>, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    for result in reader.records() {
        let record = result.context("read FASTA record")?;
        let name = record.name().to_vec();
        let seq = record.sequence().as_ref().to_vec();
        out.insert(name, seq);
    }
    Ok(out)
}

fn compute_md_nm(rec: &RecordBuf, refseq: &[u8], ref_start_1: usize) -> Option<(String, u32)> {
    use std::fmt::Write;
    let seq = rec.sequence().as_ref();
    if seq.is_empty() {
        return None;
    }
    let mut md = String::new();
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
    Some((md, nm))
}
