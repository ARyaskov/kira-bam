//! Region-restricted random access over an indexed coordinate-sorted BAM.
//!
//! `mpileup -r` normally streams the whole file and filters positions through
//! the region cursor. When a valid `.bai`/`.csi` exists, this source instead
//! walks only the BGZF chunks of the contigs the regions name — the same
//! chunks noodles' `bam::io::Reader::query` walks, but held under our
//! ownership so a pileup sweep can pull one record at a time from the merge of
//! several inputs. No record is ever buffered ahead of the pileup's sliding
//! window, whether a contig is deep or not.
//!
//! The query is per contig, not per region interval: under `-a` the sweep's
//! zero-fill bookkeeping depends on every record a contig holds (coverage the
//! regions don't admit still positions heads and contig switches), so
//! whole-interval chunk lists plus re-read deduplication could not reproduce
//! the streaming sweep exactly. Contigs no region names are skipped entirely,
//! which is the one unambiguous win — a `-r` run never decompresses the
//! untargeted majority of the file, and the fetched record stream for the
//! targeted contigs is exactly what streaming them yields.

use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_bgzf::{self as bgzf, io::Seek as _};
use noodles_core::{Position, region::Interval};
use noodles_csi::{
    self as csi, BinningIndex, binning_index::index::reference_sequence::bin::Chunk,
};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

use crate::depth::CursorMap;
use crate::io::{InputFormat, RecordSource, detect_format};

pub(crate) struct IndexedInput {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
    index: Box<dyn BinningIndex>,
    header: Header,
    /// `(local tid, whole-contig span)` for each contig `-r` names, in the
    /// input's own @SQ order — the order a coordinate-sorted file is sorted in.
    plan: Vec<(usize, Interval)>,
    plan_at: usize,
    /// Chunk queue of the contig currently being fetched; `chunk_end` is
    /// `Some` while the reader stands inside a chunk.
    chunks: std::vec::IntoIter<Chunk>,
    chunk_end: Option<bgzf::VirtualPosition>,
    path: std::path::PathBuf,
}

impl IndexedInput {
    /// Open `path` as a random-access source restricted to the contigs
    /// `regions` names, or `Ok(None)` when the input is not region-queryable
    /// (no index, or not a BGZF BAM) — such inputs keep the streaming reader.
    /// A present-but-unparsable index is an error: silently falling back to a
    /// full scan would run the whole sweep believing records are missing.
    ///
    /// `regions` is the `-r` cursor map `mpileup` already parsed against the
    /// union header; only its contig names drive eligibility here.
    pub(crate) fn open(path: &Path, regions: &CursorMap) -> Result<Option<Self>> {
        if detect_format(path).with_context(|| format!("detect format {}", path.display()))?
            != InputFormat::Bam
        {
            return Ok(None);
        }
        let Some(index) = read_associated_index(path)? else {
            return Ok(None);
        };

        let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let mut reader = bam::io::Reader::from(bgzf::io::Reader::new(file));
        let header = reader.read_header().context("read BAM header")?;

        let mut plan = Vec::new();
        for (tid, (name, map)) in header.reference_sequences().iter().enumerate() {
            if !regions.contains_key(&name.to_string()) {
                continue;
            }
            plan.push((tid, interval(1, usize::from(map.length()) as u64)?));
        }

        Ok(Some(Self {
            reader,
            index,
            header,
            plan,
            plan_at: 0,
            chunks: std::vec::Vec::<Chunk>::new().into_iter(),
            chunk_end: None,
            path: path.to_path_buf(),
        }))
    }
}

/// `.[b|c]si` path next to a BAM, as noodles' `indexed_reader::Builder` and
/// `kira-bam index` both spell it.
fn with_ext(path: &Path, ext: &str) -> std::path::PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    std::path::PathBuf::from(s)
}

/// `.bai` first, then `.csi` — the resolution order noodles uses. `None` when
/// neither index file exists; a readable but invalid one is an error with an
/// actionable rebuild hint, never a silent fallback.
fn read_associated_index(path: &Path) -> Result<Option<Box<dyn BinningIndex>>> {
    let bai_path = with_ext(path, "bai");
    if bai_path.exists() {
        let index: bam::bai::Index = bam::bai::fs::read(&bai_path).with_context(|| {
            format!(
                "index {} is corrupt or stale — delete it and rebuild with `kira-bam index {}`",
                bai_path.display(),
                path.display()
            )
        })?;
        return Ok(Some(Box::new(index)));
    }
    let csi_path = with_ext(path, "csi");
    if csi_path.exists() {
        let index: csi::Index = csi::fs::read(&csi_path).with_context(|| {
            format!(
                "index {} is corrupt or stale — delete it and rebuild with `kira-bam index {}`",
                csi_path.display(),
                path.display()
            )
        })?;
        return Ok(Some(Box::new(index)));
    }
    Ok(None)
}

/// 1-based inclusive `[start, end]` reference span as a noodles interval.
fn interval(start: u64, end: u64) -> Result<Interval> {
    let s = Position::try_from(start as usize)
        .map_err(|e| anyhow::anyhow!("region start {start} out of range: {e}"))?;
    let e = Position::try_from(end as usize)
        .map_err(|e| anyhow::anyhow!("region end {end} out of range: {e}"))?;
    Ok((s..=e).into())
}

impl RecordSource for IndexedInput {
    fn read_record_buf(&mut self, out: &mut RecordBuf) -> Result<bool> {
        loop {
            if let Some(end) = self.chunk_end {
                if self.reader.get_ref().virtual_position() >= end {
                    self.chunk_end = None;
                    continue;
                }
                let n = self
                    .reader
                    .read_record_buf(&self.header, out)
                    .with_context(|| format!("read record from {}", self.path.display()))?;
                if n == 0 {
                    self.chunk_end = None;
                    continue;
                }
                return Ok(true);
            }
            if let Some(chunk) = self.chunks.next() {
                self.reader
                    .get_mut()
                    .seek_to_virtual_position(chunk.start())
                    .with_context(|| {
                        format!(
                            "seek {}, the index may be stale — delete it and rebuild with \
                             `kira-bam index {}`",
                            self.path.display(),
                            self.path.display()
                        )
                    })?;
                self.chunk_end = Some(chunk.end());
                continue;
            }
            // Next planned contig; `false` when every named contig is done.
            if !self.advance()? {
                return Ok(false);
            }
        }
    }
}

impl IndexedInput {
    /// Query the index for the next planned contig's chunks. `false` when the
    /// plan is exhausted.
    fn advance(&mut self) -> Result<bool> {
        while self.plan_at < self.plan.len() {
            let (tid, iv) = &self.plan[self.plan_at];
            self.plan_at += 1;
            let chunks = self
                .index
                .query(*tid, iv.clone())
                .with_context(|| format!("query index of {}", self.path.display()))?;
            if !chunks.is_empty() {
                self.chunks = chunks.into_iter();
                return Ok(true);
            }
        }
        Ok(false)
    }
}
