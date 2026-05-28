use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;

use crate::cli::CatArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg};
use crate::types::OutputFormat;

pub fn run(args: CatArgs) -> Result<()> {
    let first_input = args.header_from.as_ref().unwrap_or_else(|| &args.inputs[0]);
    let header_reader = BamReader::open(first_input).context("open header source")?;
    let mut header = header_reader.header().clone();
    drop(header_reader);
    append_pg(&mut header, &PgInfo::new("cat", !args.no_pg))?;

    let fmt = if args.compression == 0 {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    };
    let mut writer = BamWriter::create(Some(&args.output), header.clone(), fmt)?;
    writer.write_header()?;

    let mut rec = RecordBuf::default();
    for input in &args.inputs {
        let mut reader =
            BamReader::open(input).with_context(|| format!("open {}", input.display()))?;
        let _ = reader.header();
        while reader.read_record_buf(&mut rec)? {
            let _ = writer.write_record_buf(&rec);
        }
    }
    writer.finish()?;
    Ok(())
}
