use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::RecordBuf;
use rustc_hash::FxHasher;

use crate::cli::SplitArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

pub fn run(args: SplitArgs) -> Result<()> {
    install_thread_pool(args.threads);

    if args.tag.len() != 2 {
        anyhow::bail!("tag must be 2 chars");
    }
    let tag_bytes = args.tag.as_bytes();
    let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

    let mut reader = BamReader::open(&args.input).context("open input")?;
    let header = reader.header().clone();
    let input_stem = args
        .input
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("split");

    let mut writers: HashMap<String, BamWriter, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();

    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        let key = extract_tag_string(&rec, tag).unwrap_or_else(|| b"unknown".to_vec());
        let key_str = String::from_utf8_lossy(&key).into_owned();
        let writer = if let Some(w) = writers.get_mut(&key_str) {
            w
        } else {
            let out = expand_format(&args.format, input_stem, &key_str);
            let mut hh = header.clone();
            append_pg(&mut hh, &PgInfo::new("split", !args.no_pg))?;
            let mut w = BamWriter::create(Some::<PathBuf>(out.into()), hh, OutputFormat::Bam)?;
            w.write_header()?;
            writers.insert(key_str.clone(), w);
            writers.get_mut(&key_str).unwrap()
        };
        let _ = writer.write_record_buf(&rec);
    }
    for (_, w) in writers {
        w.finish()?;
    }
    Ok(())
}

fn extract_tag_string(rec: &RecordBuf, t: Tag) -> Option<Vec<u8>> {
    match rec.data().get(&t) {
        Some(Value::String(s)) => Some(s.to_vec()),
        Some(Value::Character(c)) => Some(vec![*c]),
        _ => None,
    }
}

fn expand_format(fmt: &str, stem: &str, tag_value: &str) -> String {
    fmt.replace("%*", stem).replace("%!", tag_value)
}
