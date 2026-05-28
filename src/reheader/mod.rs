use std::fs::File;

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

use crate::cli::ReheaderArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg};
use crate::types::OutputFormat;

pub fn run(args: ReheaderArgs) -> Result<()> {
    let mut hdr_file = File::open(&args.header).context("open new header file")?;
    let mut sam_reader = noodles_sam::io::Reader::new(std::io::BufReader::new(&mut hdr_file));
    let mut new_header: Header = sam_reader.read_header().context("parse SAM header")?;
    append_pg(&mut new_header, &PgInfo::new("reheader", !args.no_pg))?;

    let mut reader = BamReader::open(&args.input).context("open input")?;
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
    while reader.read_record_buf(&mut rec)? {
        let _ = writer.write_record_buf(&rec);
    }
    writer.finish()?;
    Ok(())
}
