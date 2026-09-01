use std::fs::File;
use std::io::BufWriter;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_csi as csi;
use noodles_csi::binning_index::Indexer;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::index::{BinnedIndex, LinearIndex};
use noodles_sam::alignment::Record as _;

use crate::cli::IndexArgs;
use crate::io::default_index_path;

pub fn run(args: IndexArgs) -> Result<()> {
    let csi_format = args.csi;
    let out = args
        .output
        .clone()
        .unwrap_or_else(|| default_index_path(&args.input, csi_format));
    if csi_format {
        build_csi(&args, &out)
    } else {
        build_bai(&args, &out)
    }
}

fn build_bai(args: &IndexArgs, out: &std::path::Path) -> Result<()> {
    let mut indexer: Indexer<LinearIndex> = Indexer::default();
    let n_refs = index_records(args, &mut indexer)?;
    let index = indexer.build(n_refs);
    bam::bai::fs::write(out, &index).context("write BAI")?;
    Ok(())
}

fn build_csi(args: &IndexArgs, out: &std::path::Path) -> Result<()> {
    let mut indexer: Indexer<BinnedIndex> = Indexer::new(args.min_shift, args.depth);
    let n_refs = index_records(args, &mut indexer)?;
    let index = indexer.build(n_refs);

    let file = File::create(out).context("create .csi")?;
    let mut writer = csi::io::Writer::new(BufWriter::with_capacity(1 << 20, file));
    writer.write_index(&index).context("write CSI")?;
    Ok(())
}

/// Feed every record to `indexer`, returning the header's reference count.
///
/// Both the header's declared sort order *and* the actual record order are
/// checked. The header alone is not enough: a BAM can carry `SO:coordinate`
/// and still have positions out of order within a contig, and the index built
/// from it is wrong in a way nothing reports — queries silently miss reads.
fn index_records<I>(args: &IndexArgs, indexer: &mut Indexer<I>) -> Result<usize>
where
    Indexer<I>: RecordSink,
{
    let file = File::open(&args.input).with_context(|| format!("open {}", args.input.display()))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().context("read header")?;
    check_sort_order_header(&header)?;

    let mut record = bam::Record::default();
    let mut start_position = reader.get_ref().virtual_position();
    let mut last: Option<(usize, usize)> = None;

    while reader.read_record(&mut record).context("read record")? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);
        let alignment_ctx = alignment_context(&record).context("alignment context")?;

        if let Some((tid, start, _, _)) = alignment_ctx {
            let start = usize::from(start);
            if let Some((last_tid, last_start)) = last
                && (tid < last_tid || (tid == last_tid && start < last_start))
            {
                let name = |i: usize| {
                    header
                        .reference_sequences()
                        .get_index(i)
                        .map(|(n, _)| n.to_string())
                        .unwrap_or_else(|| format!("#{i}"))
                };
                anyhow::bail!(
                    "{} is not coordinate-sorted: {}:{start} follows {}:{last_start}. \
                     Indexing it would produce a .bai/.csi that silently misses reads — \
                     sort it first (`kira-bam sort`).",
                    args.input.display(),
                    name(tid),
                    name(last_tid)
                );
            }
            last = Some((tid, start));
        }

        indexer.add(alignment_ctx, chunk)?;
        start_position = end_position;
    }
    Ok(header.reference_sequences().len())
}

/// `Indexer` is generic over its index kind but `add_record` is defined per
/// instantiation; this is just enough of a trait to share the loop above.
trait RecordSink {
    fn add(
        &mut self,
        ctx: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> std::io::Result<()>;
}

impl RecordSink for Indexer<LinearIndex> {
    fn add(
        &mut self,
        ctx: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> std::io::Result<()> {
        self.add_record(ctx, chunk)
    }
}

impl RecordSink for Indexer<BinnedIndex> {
    fn add(
        &mut self,
        ctx: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> std::io::Result<()> {
        self.add_record(ctx, chunk)
    }
}

fn check_sort_order_header(header: &noodles_sam::Header) -> Result<()> {
    use noodles_sam::header::record::value::map::header::tag::SORT_ORDER;

    let so = header
        .header()
        .and_then(|hd| hd.other_fields().get(&SORT_ORDER))
        .map(|v| v.to_vec())
        .unwrap_or_default();
    if so != b"coordinate" {
        let shown = if so.is_empty() {
            "unspecified".to_string()
        } else {
            String::from_utf8_lossy(&so).into_owned()
        };
        anyhow::bail!(
            "@HD sort order is {shown}, not coordinate. Only a coordinate-sorted BAM can be \
             indexed — sort it first (`kira-bam sort`)."
        );
    }
    Ok(())
}

fn alignment_context(
    record: &bam::Record,
) -> std::io::Result<Option<(usize, Position, Position, bool)>> {
    let tid = record.reference_sequence_id().transpose()?;
    let start = record.alignment_start().transpose()?;
    let end = record.alignment_end().transpose()?;
    match (tid, start, end) {
        (Some(id), Some(s), Some(e)) => Ok(Some((id, s, e, !record.flags().is_unmapped()))),
        _ => Ok(None),
    }
}
