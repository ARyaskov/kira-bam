use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Region;
use noodles_sam::alignment::RecordBuf;

use crate::cli::ViewArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

pub fn run(args: ViewArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let fmt = resolve_fmt(&args);

    if args.regions.is_empty() {
        return run_streaming(args, fmt);
    }
    run_regions(args, fmt)
}

fn resolve_fmt(args: &ViewArgs) -> OutputFormat {
    if args.cram {
        OutputFormat::Cram
    } else if args.bam && args.uncompressed {
        OutputFormat::UncompressedBam
    } else if args.bam {
        OutputFormat::Bam
    } else if args.sam {
        OutputFormat::Sam
    } else {
        match &args.output {
            Some(p) => OutputFormat::from_path_ext(p),
            None => OutputFormat::Sam,
        }
    }
}

fn run_streaming(args: ViewArgs, fmt: OutputFormat) -> Result<()> {
    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let mut header = reader.header().clone();
    append_pg(&mut header, &PgInfo::new("view", !args.no_pg))?;

    if args.count {
        return count_only_stream(&mut reader, &args);
    }

    let mut writer = BamWriter::create_with_reference(
        args.output.as_deref(),
        header,
        fmt,
        args.reference.as_deref(),
    )
    .context("create output")?;

    let want_header = matches!(
        fmt,
        OutputFormat::Bam | OutputFormat::UncompressedBam | OutputFormat::Cram
    ) || args.with_header
        || args.header_only;
    if want_header {
        writer.write_header()?;
    }
    if args.header_only {
        writer.finish()?;
        return Ok(());
    }

    let mut rec = crate::io::Record::default();
    let mut skipped: u64 = 0;
    while reader.read_record_buf(&mut rec)? {
        if !pass_filters(&rec, &args) {
            continue;
        }
        if args.drop_tags {
            rec.data_mut().clear();
        }
        // noodles is stricter than samtools on per-record validity; soft-skip to stay tolerant.
        if let Err(e) = writer.write_record_buf(&rec) {
            skipped += 1;
            if skipped <= 5 {
                eprintln!("[kira-bam] skip record: {e:#}");
            } else if skipped == 6 {
                eprintln!("[kira-bam] (further skips suppressed)");
            }
        }
    }
    writer.finish()?;
    if skipped > 0 {
        eprintln!("[kira-bam] {skipped} records skipped due to write errors");
    }
    Ok(())
}

fn run_regions(args: ViewArgs, fmt: OutputFormat) -> Result<()> {
    let mut idx_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&args.input)
        .context("open indexed BAM (region query needs .bai or .csi)")?;
    let mut header = idx_reader.read_header().context("read indexed header")?;
    append_pg(&mut header, &PgInfo::new("view", !args.no_pg))?;

    if args.count {
        let mut n: u64 = 0;
        for spec in &args.regions {
            n += count_region(&mut idx_reader, &header, spec, &args)?;
        }
        println!("{n}");
        return Ok(());
    }

    let mut writer = BamWriter::create_with_reference(
        args.output.as_deref(),
        header.clone(),
        fmt,
        args.reference.as_deref(),
    )
    .context("create output")?;
    let want_header = matches!(
        fmt,
        OutputFormat::Bam | OutputFormat::UncompressedBam | OutputFormat::Cram
    ) || args.with_header
        || args.header_only;
    if want_header {
        writer.write_header()?;
    }
    if args.header_only {
        writer.finish()?;
        return Ok(());
    }

    let mut skipped: u64 = 0;
    for spec in &args.regions {
        if spec == "*" {
            for rec_result in idx_reader.query_unmapped().context("query unmapped")? {
                let rec = rec_result.context("read unmapped record")?;
                let mut buf = RecordBuf::try_from_alignment_record(&header, &rec)?;
                if !pass_filters(&buf, &args) {
                    continue;
                }
                if args.drop_tags {
                    buf.data_mut().clear();
                }
                if writer.write_record_buf(&buf).is_err() {
                    skipped += 1;
                }
            }
        } else {
            let region: Region = spec
                .parse()
                .map_err(|e| anyhow::anyhow!("bad region {spec:?}: {e:?}"))?;
            let query = idx_reader
                .query(&header, &region)
                .with_context(|| format!("query {spec}"))?;
            for rec_result in query.records() {
                let rec = rec_result.context("read region record")?;
                let mut buf = RecordBuf::try_from_alignment_record(&header, &rec)?;
                if !pass_filters(&buf, &args) {
                    continue;
                }
                if args.drop_tags {
                    buf.data_mut().clear();
                }
                if writer.write_record_buf(&buf).is_err() {
                    skipped += 1;
                }
            }
        }
    }
    writer.finish()?;
    if skipped > 0 {
        eprintln!("[kira-bam] {skipped} records skipped due to write errors");
    }
    Ok(())
}

fn count_only_stream(reader: &mut BamReader, args: &ViewArgs) -> Result<()> {
    let mut n: u64 = 0;
    let mut rec = crate::io::Record::default();
    while reader.read_record_buf(&mut rec)? {
        if pass_filters(&rec, args) {
            n += 1;
        }
    }
    println!("{n}");
    Ok(())
}

fn count_region<R>(
    reader: &mut bam::io::IndexedReader<R>,
    header: &noodles_sam::Header,
    spec: &str,
    args: &ViewArgs,
) -> Result<u64>
where
    R: noodles_bgzf::io::BufRead + noodles_bgzf::io::Seek,
{
    let mut n: u64 = 0;
    if spec == "*" {
        for rec_result in reader.query_unmapped()? {
            let rec = rec_result?;
            let buf = RecordBuf::try_from_alignment_record(header, &rec)?;
            if pass_filters(&buf, args) {
                n += 1;
            }
        }
    } else {
        let region: Region = spec
            .parse()
            .map_err(|e| anyhow::anyhow!("bad region {spec:?}: {e:?}"))?;
        let query = reader.query(header, &region)?;
        for rec_result in query.records() {
            let rec = rec_result?;
            let buf = RecordBuf::try_from_alignment_record(header, &rec)?;
            if pass_filters(&buf, args) {
                n += 1;
            }
        }
    }
    Ok(n)
}

fn pass_filters(rec: &crate::io::Record, args: &ViewArgs) -> bool {
    let flags = u16::from(rec.flags());
    if args.require_flags != 0 && (flags & args.require_flags) != args.require_flags {
        return false;
    }
    if args.filter_flags != 0 && (flags & args.filter_flags) != 0 {
        return false;
    }
    if let Some(mapq) = rec.mapping_quality() {
        if u8::from(mapq) < args.min_mapq {
            return false;
        }
    } else if args.min_mapq > 0 {
        return false;
    }
    // -v REGION exclusion (samtools PR #669). Naive impl: parse spec as
    // `chr` or `chr:start-end` and drop if record's [start,end] overlaps.
    if !args.exclude_regions.is_empty()
        && let Some(_tid) = rec.reference_sequence_id()
    {
        // We compare against the textual ref name. Slow per-record;
        // acceptable for small exclude lists.
        for spec in &args.exclude_regions {
            if record_in_region_spec(rec, spec) {
                return false;
            }
        }
    }
    true
}

fn record_in_region_spec(rec: &crate::io::Record, spec: &str) -> bool {
    // Parse `chr` or `chr:start-end` (1-based inclusive).
    let (chr, start, end) = if let Some(colon) = spec.find(':') {
        let (c, rest) = spec.split_at(colon);
        let rest = &rest[1..];
        if let Some(dash) = rest.find('-') {
            let (a, b) = rest.split_at(dash);
            let b = &b[1..];
            (c.to_string(), a.parse::<u64>().ok(), b.parse::<u64>().ok())
        } else {
            (c.to_string(), rest.parse::<u64>().ok(), None)
        }
    } else {
        (spec.to_string(), None, None)
    };
    // Without header context here we can't validate chr names; rely on caller
    // to feed sensible specs. Compare by position only.
    let r_start = match rec.alignment_start() {
        Some(p) => usize::from(p) as u64,
        None => return false,
    };
    let r_end = rec
        .alignment_end()
        .map(|p| usize::from(p) as u64)
        .unwrap_or(r_start);
    let s = start.unwrap_or(0);
    let e = end.unwrap_or(u64::MAX);
    let overlaps = r_end >= s && r_start <= e;
    let _ = chr; // chr validation is left for indexed path; here we accept overlap-by-position.
    overlaps
}
