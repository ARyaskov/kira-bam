use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use rustc_hash::FxHasher;

use crate::cli::SplitArgs;
use crate::io::{
    BamReader, BamWriter, PgInfo, append_pg, ensure_distinct_output, report_skipped,
    write_record_tolerant,
};
use crate::types::OutputFormat;

#[derive(Debug, Eq, Hash, PartialEq)]
enum SplitKey {
    Missing,
    Text(Vec<u8>),
    Integer(i64),
    Float(u32),
    Hex(Vec<u8>),
}

impl SplitKey {
    fn from_value(value: Option<&Value>) -> Result<Self> {
        Ok(match value {
            None => Self::Missing,
            Some(Value::String(s)) => Self::Text(s.to_vec()),
            Some(Value::Character(c)) => Self::Text(vec![*c]),
            Some(Value::Float(n)) => Self::Float(n.to_bits()),
            Some(Value::Hex(s)) => Self::Hex(s.to_vec()),
            Some(Value::Array(_)) => anyhow::bail!("array-valued tags cannot be split"),
            Some(value) => Self::Integer(value.as_int().context("expected an integer tag value")?),
        })
    }

    fn filename_component(&self) -> String {
        // Reserve 'unknown' for missing tags. Type prefixes cannot be
        // produced by escape_component, which encodes '%' as '%25'.
        match self {
            Self::Missing => "unknown".to_owned(),
            Self::Text(s) if s == b"unknown" => "%75nknown".to_owned(),
            Self::Text(s) => escape_component(s),
            Self::Integer(n) => format!("%i{n}"),
            Self::Float(bits) => format!("%f_{bits:08X}"),
            Self::Hex(s) => format!("%h{}", escape_component(s)),
        }
    }
}

pub fn run(args: SplitArgs) -> Result<()> {
    if args.tag.len() != 2 {
        anyhow::bail!("tag must be 2 chars");
    }
    let tag_bytes = args.tag.as_bytes();
    let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let header = reader.header().clone();
    let input_stem = args
        .input
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("split");

    let mut writers: HashMap<SplitKey, BamWriter, std::hash::BuildHasherDefault<FxHasher>> =
        HashMap::default();
    let mut destinations = HashSet::new();

    let mut rec = RecordBuf::default();
    let mut skipped = 0u64;
    while reader.read_record_buf(&mut rec)? {
        let key = SplitKey::from_value(rec.data().get(&tag))
            .with_context(|| format!("split on tag {}", args.tag))?;
        let writer = if let Some(w) = writers.get_mut(&key) {
            w
        } else {
            // A BAM writer cannot be reopened once its BGZF EOF block is out, so
            // every distinct value has to stay open for the whole pass. Refuse
            // loudly at the cap rather than let the OS report a bare
            // "too many open files" halfway through.
            if writers.len() >= args.max_outputs {
                anyhow::bail!(
                    "tag {} has more than {} distinct values (at {key:?}); each needs its own \
                     open output file. Raise --max-outputs if your fd limit allows, or split \
                     on a lower-cardinality tag.",
                    args.tag,
                    args.max_outputs
                );
            }
            let out = PathBuf::from(expand_format(
                &args.format,
                input_stem,
                &key.filename_component(),
            ));
            anyhow::ensure!(
                out != std::path::Path::new("-"),
                "split expects file destinations"
            );
            let destination = ensure_distinct_output(&args.input, &out)?;
            anyhow::ensure!(
                destinations.insert(destination),
                "multiple tag values resolve to {}",
                out.display()
            );
            let mut hh = header.clone();
            append_pg(&mut hh, &PgInfo::new("split", !args.no_pg))?;
            let mut w = BamWriter::create(Some(out), hh, OutputFormat::Bam)?;
            w.write_header()?;
            writers.entry(key).or_insert(w)
        };
        write_record_tolerant(writer, &rec, "split", &mut skipped)?;
    }
    for (_, w) in writers {
        w.finish()?;
    }
    report_skipped("split", skipped);
    Ok(())
}

fn expand_format(fmt: &str, stem: &str, tag_value: &str) -> String {
    let mut result = String::with_capacity(fmt.len());
    let mut characters = fmt.chars().peekable();
    while let Some(c) = characters.next() {
        match (c, characters.peek()) {
            ('%', Some('*')) => {
                characters.next();
                result.push_str(stem);
            }
            ('%', Some('!')) => {
                characters.next();
                result.push_str(tag_value);
            }
            _ => result.push(c),
        }
    }
    result
}

fn escape_component(value: &[u8]) -> String {
    use std::fmt::Write as _;
    let mut encoded = String::new();
    for &byte in value {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-') {
            encoded.push(char::from(byte));
        } else {
            write!(encoded, "%{byte:02X}").expect("writing to a String cannot fail");
        }
    }
    encoded
}
