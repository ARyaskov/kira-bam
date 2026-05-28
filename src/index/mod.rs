use std::fs::File;
use std::io::BufWriter;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_csi as csi;
use noodles_csi::binning_index::Indexer;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::index::BinnedIndex;
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
    let index = bam::fs::index(&args.input).context("build BAI index")?;
    bam::bai::fs::write(out, &index).context("write BAI")?;
    Ok(())
}

fn build_csi(args: &IndexArgs, out: &std::path::Path) -> Result<()> {
    let file = File::open(&args.input).context("open BAM")?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().context("read header")?;
    let mut indexer: Indexer<BinnedIndex> = Indexer::new(args.min_shift, args.depth);

    let mut record = bam::Record::default();
    let mut start_position = reader.get_ref().virtual_position();
    while reader.read_record(&mut record).context("read record")? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);
        let alignment_ctx = alignment_context(&record).context("alignment context")?;
        indexer.add_record(alignment_ctx, chunk)?;
        start_position = end_position;
    }
    let index = indexer.build(header.reference_sequences().len());

    let file = File::create(out).context("create .csi")?;
    let mut writer = csi::io::Writer::new(BufWriter::with_capacity(1 << 20, file));
    writer.write_index(&index).context("write CSI")?;
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
