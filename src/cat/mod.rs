//! Concatenate BAMs by copying their compressed blocks.
//!
//! Only the output header, and the tail of an input's header block when it does
//! not end on a block boundary, are compressed; the rest is copied byte for byte.

use std::fs::File;
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::Header;

use crate::cli::CatArgs;
use crate::io::{BGZF_EOF, BamReader, PgInfo, append_pg};

pub fn run(args: CatArgs) -> Result<()> {
    let header_src = args.header_from.as_ref().unwrap_or(&args.inputs[0]);
    let reader = BamReader::open_with_reference(header_src, args.reference.as_deref())
        .context("open header source")?;
    let mut header = reader.header().clone();
    drop(reader);
    append_pg(&mut header, &PgInfo::new("cat", !args.no_pg))?;

    // Checked up front so a mismatch leaves no half-written output.
    for input in &args.inputs {
        let r = BamReader::open_with_reference(input, args.reference.as_deref())
            .with_context(|| format!("open {}", input.display()))?;
        check_same_references(&header, r.header(), input)?;
    }

    let file = File::create(&args.output).context("create output")?;
    let mut sink = BufWriter::with_capacity(1 << 20, file);

    let level = compression_level(args.compression);
    write_header_block(&mut sink, &header, level).context("write output header")?;

    for input in &args.inputs {
        copy_records(input, &mut sink, level)
            .with_context(|| format!("copy records from {}", input.display()))?;
    }

    sink.write_all(&BGZF_EOF).context("write BGZF EOF block")?;
    sink.flush().context("flush output")?;
    Ok(())
}

fn compression_level(level: i32) -> Option<bgzf::io::writer::CompressionLevel> {
    if level < 0 {
        None
    } else {
        bgzf::io::writer::CompressionLevel::new(level.clamp(0, 9) as u8)
    }
}

/// Write the header into whole BGZF blocks, leaving the sink at a block boundary.
fn write_header_block<W: Write>(
    sink: &mut W,
    header: &Header,
    level: Option<bgzf::io::writer::CompressionLevel>,
) -> Result<()> {
    let mut bgz = new_bgzf_writer(sink, level);
    let mut w = bam::io::Writer::from(&mut bgz);
    w.write_header(header)?;
    // `flush` closes the block; `try_finish` would also append an EOF marker,
    // which mid-stream ends the file early.
    bgz.flush()?;
    Ok(())
}

fn new_bgzf_writer<W: Write>(
    sink: W,
    level: Option<bgzf::io::writer::CompressionLevel>,
) -> bgzf::io::Writer<W> {
    let mut builder = bgzf::io::writer::Builder::default();
    if let Some(level) = level {
        builder = builder.set_compression_level(level);
    }
    builder.build_from_writer(sink)
}

/// Append one input's record stream to `sink`.
fn copy_records<W: Write>(
    path: &Path,
    sink: &mut W,
    level: Option<bgzf::io::writer::CompressionLevel>,
) -> Result<()> {
    use std::io::BufRead as _;

    let file = File::open(path)?;
    let bgz = bgzf::io::Reader::new(file);
    let mut reader = bam::io::Reader::from(bgz);
    reader.read_header().context("read header")?;
    let mut bgz = reader.into_inner();

    // A header rarely ends on a block boundary, and the record bytes sharing its
    // last block are the one part that cannot be copied.
    let mut resume = bgz.virtual_position();
    if resume.uncompressed() > 0 {
        let tail = bgz.fill_buf()?.to_vec();
        let after = bgz.virtual_position();
        if after.compressed() == resume.compressed() {
            bgz.consume(tail.len());
            if !tail.is_empty() {
                let mut out = new_bgzf_writer(&mut *sink, level);
                out.write_all(&tail)?;
                out.flush()?;
            }
            // Loads the next block, whose position is where copying resumes.
            bgz.fill_buf()?;
            resume = bgz.virtual_position();
        } else {
            // Already exhausted: `fill_buf` moved on and we are block-aligned.
            resume = after;
        }
    }
    let start = resume.compressed();
    drop(bgz);

    copy_raw_range(path, start, sink)?;
    Ok(())
}

/// Copy `path` from `start` to just before its BGZF EOF marker, verbatim.
fn copy_raw_range<W: Write>(path: &Path, start: u64, sink: &mut W) -> Result<u64> {
    let mut file = File::open(path)?;
    let len = file.metadata()?.len();

    // Record blocks are exactly the bytes before the trailing marker.
    let mut tail = [0u8; BGZF_EOF.len()];
    if len < tail.len() as u64 {
        anyhow::bail!("shorter than a BGZF EOF marker — truncated?");
    }
    file.seek(SeekFrom::End(-(tail.len() as i64)))?;
    file.read_exact(&mut tail)?;
    if tail != BGZF_EOF {
        anyhow::bail!("missing BGZF EOF marker — truncated or not BGZF");
    }
    let end = len - tail.len() as u64;
    if end <= start {
        return Ok(0);
    }

    file.seek(SeekFrom::Start(start))?;
    let mut src = file.take(end - start);
    io::copy(&mut src, sink).context("copy compressed blocks")
}

/// Reference ids are copied untouched, so a differing @SQ list would relabel
/// every alignment.
fn check_same_references(out: &Header, input: &Header, path: &Path) -> Result<()> {
    let a = out.reference_sequences();
    let b = input.reference_sequences();
    let same = a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|((an, aq), (bn, bq))| an == bn && aq.length() == bq.length());
    if same {
        return Ok(());
    }
    anyhow::bail!(
        "{}: @SQ list differs from the output header, so its reference ids would be \
         reinterpreted. Use `kira-bam merge` instead, or reheader the inputs to match.",
        path.display()
    )
}
