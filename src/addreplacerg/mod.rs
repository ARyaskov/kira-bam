use anyhow::{Context, Result};
use bstr::BString;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::header::record::value::{Map, map::ReadGroup, map::read_group::tag as rg_tag};

use crate::cli::AddReplaceRgArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

pub fn run(args: AddReplaceRgArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let mut reader = BamReader::open(&args.input).context("open input")?;
    let mut header = reader.header().clone();
    let (rg_id, rg_map) = parse_rg_line(&args.rg_line)?;
    header
        .read_groups_mut()
        .insert(rg_id.clone(), rg_map);
    append_pg(&mut header, &PgInfo::new("addreplacerg", !args.no_pg))?;

    let fmt = if args.uncompressed {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    };
    let mut writer = BamWriter::create(Some(&args.output), header, fmt)?;
    writer.write_header()?;

    let mode = args.mode.as_str();
    let mut rec = RecordBuf::default();
    while reader.read_record_buf(&mut rec)? {
        let has_rg = rec.data().get(&Tag::READ_GROUP).is_some();
        let should_set = match mode {
            "overwrite_all" => true,
            "overwrite_orphans" | "orphan-only" => !has_rg,
            "lenient" => !has_rg,
            _ => !has_rg,
        };
        if should_set {
            rec.data_mut()
                .insert(Tag::READ_GROUP, Value::String(rg_id.clone()));
        }
        let _ = writer.write_record_buf(&rec);
    }
    writer.finish()?;
    Ok(())
}

fn parse_rg_line(line: &str) -> Result<(BString, Map<ReadGroup>)> {
    // Accept either `@RG\tID:foo\tSM:bar` or `ID:foo\tSM:bar`.
    let body = line.strip_prefix("@RG\t").unwrap_or(line);
    let mut id: Option<BString> = None;
    let mut builder = Map::<ReadGroup>::builder();
    for field in body.split('\t') {
        let (k, v) = field.split_once(':').context("malformed @RG field")?;
        match k {
            "ID" => id = Some(BString::from(v.as_bytes())),
            "SM" => builder = builder.insert(rg_tag::SAMPLE, v.as_bytes()),
            "LB" => builder = builder.insert(rg_tag::LIBRARY, v.as_bytes()),
            "PL" => builder = builder.insert(rg_tag::PLATFORM, v.as_bytes()),
            "PU" => builder = builder.insert(rg_tag::PLATFORM_UNIT, v.as_bytes()),
            "CN" => builder = builder.insert(rg_tag::SEQUENCING_CENTER, v.as_bytes()),
            "DT" => builder = builder.insert(rg_tag::PRODUCED_AT, v.as_bytes()),
            "PI" => builder = builder.insert(rg_tag::PREDICTED_MEDIAN_INSERT_SIZE, v.as_bytes()),
            "DS" => builder = builder.insert(rg_tag::DESCRIPTION, v.as_bytes()),
            _other => {
                // Unknown RG fields ignored; samtools tolerates but we drop for simplicity.
            }
        }
    }
    let id = id.context("missing ID in @RG line")?;
    let map = builder.build().context("build @RG")?;
    Ok((id, map))
}
