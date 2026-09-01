use std::fs::File;

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

use crate::cli::ReheaderArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, report_skipped, write_record_tolerant};
use crate::types::OutputFormat;

pub fn run(args: ReheaderArgs) -> Result<()> {
    let mut hdr_file = File::open(&args.header).context("open new header file")?;
    let mut sam_reader = noodles_sam::io::Reader::new(std::io::BufReader::new(&mut hdr_file));
    let mut new_header: Header = sam_reader.read_header().context("parse SAM header")?;
    append_pg(&mut new_header, &PgInfo::new("reheader", !args.no_pg))?;

    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    // Validate @SQ compatibility before replacing.
    let old_header = reader.header().clone();
    if old_header.reference_sequences().len() != new_header.reference_sequences().len() {
        eprintln!(
            "[kira-bam reheader] warning: @SQ count differs (old {}, new {}). Records still reference old IDs.",
            old_header.reference_sequences().len(),
            new_header.reference_sequences().len()
        );
    }

    let mut writer = BamWriter::create(Some(&args.output), new_header, OutputFormat::Bam)?;
    writer.write_header()?;
    let mut rec = RecordBuf::default();
    let mut skipped = 0u64;
    while reader.read_record_buf(&mut rec)? {
        write_record_tolerant(&mut writer, &rec, "reheader", &mut skipped)?;
    }
    writer.finish()?;
    report_skipped("reheader", skipped);
    Ok(())
}
