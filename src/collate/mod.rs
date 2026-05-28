use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;

use crate::cli::CollateArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;
use rustc_hash::FxHasher;
use std::collections::HashMap;

pub fn run(args: CollateArgs) -> Result<()> {
    install_thread_pool(args.threads);

    // Simple in-memory implementation: group by qname hash.
    // For huge files we'd spill — but this fits the README's "early development" disclaimer.

    let mut reader = BamReader::open(&args.input).context("open input")?;
    let mut header = reader.header().clone();
    append_pg(&mut header, &PgInfo::new("collate", !args.no_pg))?;

    let mut buckets: HashMap<Vec<u8>, Vec<RecordBuf>, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    let mut rec = RecordBuf::default();
    let mut total = 0u64;
    while reader.read_record_buf(&mut rec)? {
        if args.fast {
            let f = u16::from(rec.flags());
            if (f & (0x100 | 0x800)) != 0 {
                continue;
            }
        }
        let qname: Vec<u8> = rec.name().map(|n| n.to_vec()).unwrap_or_default();
        buckets.entry(qname).or_default().push(rec.clone());
        total += 1;
    }

    let mut out_path = args.output.clone();
    if out_path.extension().is_none() {
        out_path.set_extension("bam");
    }
    let mut writer = BamWriter::create(Some(&out_path), header, OutputFormat::Bam)?;
    writer.write_header()?;
    for (_q, recs) in buckets.drain() {
        for r in recs {
            let _ = writer.write_record_buf(&r);
        }
    }
    writer.finish()?;
    eprintln!(
        "[kira-bam collate] wrote {total} records to {}",
        out_path.display()
    );
    Ok(())
}
