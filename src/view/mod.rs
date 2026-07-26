use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Region;
use noodles_sam::alignment::RecordBuf;

use crate::cli::ViewArgs;
use crate::io::{
    BamReader, BamWriter, PgInfo, RecordSummary, append_pg, install_thread_pool,
    is_record_level_error,
};
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
    let exclude = resolve_exclude_regions(&args.exclude_regions, &header)?;
    append_pg(&mut header, &PgInfo::new("view", !args.no_pg))?;

    if args.count {
        return count_only_stream(&mut reader, &args, &exclude);
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
        if !pass_filters(&rec, &args, &exclude) {
            continue;
        }
        if args.drop_tags {
            rec.data_mut().clear();
        }
        // noodles is stricter than samtools on record validity; skip those, but
        // let a failing output stream through.
        if let Err(e) = writer.write_record_buf(&rec) {
            if !is_record_level_error(&e) {
                return Err(e);
            }
            skipped += 1;
            if skipped <= 5 {
                eprintln!("[kira-bam] skip record: {e:#}");
            } else if skipped == 6 {
                eprintln!("[kira-bam] (further skips suppressed)");
            }
        }
    }
    writer.finish()?;
    report_skipped(skipped);
    Ok(())
}

fn report_skipped(skipped: u64) {
    if skipped > 0 {
        eprintln!("[kira-bam] {skipped} records skipped as invalid");
    }
}

fn run_regions(args: ViewArgs, fmt: OutputFormat) -> Result<()> {
    let mut idx_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&args.input)
        .context("open indexed BAM (region query needs .bai or .csi)")?;
    let mut header = idx_reader.read_header().context("read indexed header")?;
    let exclude = resolve_exclude_regions(&args.exclude_regions, &header)?;
    append_pg(&mut header, &PgInfo::new("view", !args.no_pg))?;

    if args.count {
        let mut n: u64 = 0;
        for spec in &args.regions {
            n += count_region(&mut idx_reader, &header, spec, &args, &exclude)?;
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
    let mut emit = |buf: &mut RecordBuf, skipped: &mut u64| -> Result<()> {
        if !pass_filters(buf, &args, &exclude) {
            return Ok(());
        }
        if args.drop_tags {
            buf.data_mut().clear();
        }
        if let Err(e) = writer.write_record_buf(buf) {
            if !is_record_level_error(&e) {
                return Err(e);
            }
            *skipped += 1;
        }
        Ok(())
    };
    for spec in &args.regions {
        if spec == "*" {
            for rec_result in idx_reader.query_unmapped().context("query unmapped")? {
                let rec = rec_result.context("read unmapped record")?;
                let mut buf = RecordBuf::try_from_alignment_record(&header, &rec)?;
                emit(&mut buf, &mut skipped)?;
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
                emit(&mut buf, &mut skipped)?;
            }
        }
    }
    writer.finish()?;
    report_skipped(skipped);
    Ok(())
}

fn count_only_stream(
    reader: &mut BamReader,
    args: &ViewArgs,
    exclude: &[ExcludeRegion],
) -> Result<()> {
    let mut n: u64 = 0;
    if !exclude.is_empty() {
        // `-v` needs the alignment end, hence the CIGAR.
        let mut rec = crate::io::Record::default();
        while reader.read_record_buf(&mut rec)? {
            if pass_filters(&rec, args, exclude) {
                n += 1;
            }
        }
    } else {
        let mut rec = RecordSummary::default();
        while reader.read_summary(&mut rec)? {
            if summary_passes(&rec, args) {
                n += 1;
            }
        }
    }
    println!("{n}");
    Ok(())
}

/// The flag and MAPQ half of [`pass_filters`]. Region exclusion is excluded: it
/// needs an alignment end, which the summary does not carry.
fn summary_passes(rec: &RecordSummary, args: &ViewArgs) -> bool {
    if args.require_flags != 0 && (rec.flags & args.require_flags) != args.require_flags {
        return false;
    }
    if args.filter_flags != 0 && (rec.flags & args.filter_flags) != 0 {
        return false;
    }
    match rec.mapping_quality {
        Some(mapq) if mapq < args.min_mapq => false,
        None if args.min_mapq > 0 => false,
        _ => true,
    }
}

fn count_region<R>(
    reader: &mut bam::io::IndexedReader<R>,
    header: &noodles_sam::Header,
    spec: &str,
    args: &ViewArgs,
    exclude: &[ExcludeRegion],
) -> Result<u64>
where
    R: noodles_bgzf::io::BufRead + noodles_bgzf::io::Seek,
{
    let mut n: u64 = 0;
    if spec == "*" {
        for rec_result in reader.query_unmapped()? {
            let rec = rec_result?;
            let buf = RecordBuf::try_from_alignment_record(header, &rec)?;
            if pass_filters(&buf, args, exclude) {
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
            if pass_filters(&buf, args, exclude) {
                n += 1;
            }
        }
    }
    Ok(n)
}

/// A `-v` exclusion resolved against the header: reference id plus a 1-based
/// inclusive span.
#[derive(Clone, Copy, Debug)]
struct ExcludeRegion {
    tid: usize,
    start: u64,
    end: u64,
}

/// Parse `chr`, `chr:start` or `chr:start-end` against `header`. An unknown name
/// is an error, not a filter that silently matches nothing.
fn resolve_exclude_regions(
    specs: &[String],
    header: &noodles_sam::Header,
) -> Result<Vec<ExcludeRegion>> {
    let mut out = Vec::with_capacity(specs.len());
    for spec in specs {
        let (name, range) = match spec.rsplit_once(':') {
            Some((name, range)) => (name, Some(range)),
            None => (spec.as_str(), None),
        };
        let tid = header
            .reference_sequences()
            .get_index_of(name.as_bytes())
            .with_context(|| format!("-v {spec}: no reference sequence named {name:?}"))?;
        let (start, end) = match range {
            None => (1, u64::MAX),
            Some(range) => {
                let parse = |s: &str, what: &str| -> Result<u64> {
                    s.replace(',', "")
                        .parse::<u64>()
                        .with_context(|| format!("-v {spec}: bad {what}"))
                };
                match range.split_once('-') {
                    Some((a, b)) if !b.is_empty() => {
                        (parse(a, "start")?, parse(b, "end")?)
                    }
                    Some((a, _)) => (parse(a, "start")?, u64::MAX),
                    None => (parse(range, "start")?, u64::MAX),
                }
            }
        };
        out.push(ExcludeRegion { tid, start, end });
    }
    Ok(out)
}

fn pass_filters(rec: &crate::io::Record, args: &ViewArgs, exclude: &[ExcludeRegion]) -> bool {
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
    if !exclude.is_empty() && record_excluded(rec, exclude) {
        return false;
    }
    true
}

/// Does the record overlap any `-v` region? Unplaced records never do.
fn record_excluded(rec: &crate::io::Record, exclude: &[ExcludeRegion]) -> bool {
    let Some(tid) = rec.reference_sequence_id() else {
        return false;
    };
    let Some(start) = rec.alignment_start() else {
        return false;
    };
    let r_start = usize::from(start) as u64;
    let r_end = rec
        .alignment_end()
        .map(|p| usize::from(p) as u64)
        .unwrap_or(r_start);
    exclude
        .iter()
        .any(|r| r.tid == tid && r_end >= r.start && r_start <= r.end)
}
