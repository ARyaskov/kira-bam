//! Streaming pileup over a coordinate-sorted BAM.
//!
//! Yields `(tid, pos, [(base, qual)])` in genomic order, keeping a sliding window
//! of active reads keyed by alignment end. Nothing allocates per position or per
//! record in steady state: the column goes into a caller-owned buffer, decode
//! buffers are recycled as reads leave the window, and one record buffer is
//! reused throughout.

use anyhow::Result;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::op::Kind;

use crate::io::BamReader;

pub type PileupColumn = (usize, usize, Vec<(u8, u8)>);

#[derive(Debug, Default)]
pub struct ActiveRead {
    pub end: usize,
    /// Per-reference-position `(base, qual)` from this read's alignment start.
    pub bases: Vec<(u8, u8)>,
    /// 1-based reference start.
    pub ref_start: usize,
}

pub struct PileupIter {
    reader: BamReader,
    active: Vec<ActiveRead>,
    /// Decode buffers of reads that left the window.
    free: Vec<Vec<(u8, u8)>>,
    /// The record buffer, reused for every read.
    cur: RecordBuf,
    /// `cur` belongs to a later position or contig and must not be read over.
    held: bool,
    current_tid: Option<usize>,
    current_pos: usize,
    min_mapq: u8,
    min_baseq: u8,
}

impl PileupIter {
    pub fn new(reader: BamReader, min_mapq: u8, min_baseq: u8) -> Self {
        Self {
            reader,
            active: Vec::new(),
            free: Vec::new(),
            cur: RecordBuf::default(),
            held: false,
            current_tid: None,
            current_pos: 0,
            min_mapq,
            min_baseq,
        }
    }

    /// Make the next record available in `self.cur`. `false` at EOF.
    fn fetch_next(&mut self) -> Result<bool> {
        if self.held {
            self.held = false;
            return Ok(true);
        }
        self.reader.read_record_buf(&mut self.cur)
    }

    /// Next pileup column. [`Self::next_column_into`] avoids the allocation.
    pub fn next_column(&mut self) -> Result<Option<PileupColumn>> {
        let mut col = Vec::new();
        Ok(self
            .next_column_into(&mut col)?
            .map(|(tid, pos)| (tid, pos, col)))
    }

    /// Write the next column into `col` (cleared first), returning its position.
    pub fn next_column_into(&mut self, col: &mut Vec<(u8, u8)>) -> Result<Option<(usize, usize)>> {
        loop {
            // Pull in records that start at or before current_pos for current_tid.
            loop {
                if !self.fetch_next()? {
                    if self.active.is_empty() && self.current_tid.is_none() {
                        return Ok(None);
                    }
                    break;
                }
                let flags = u16::from(self.cur.flags());
                if (flags & 0x4) != 0 {
                    continue;
                }
                if let Some(mq) = self.cur.mapping_quality()
                    && u8::from(mq) < self.min_mapq
                {
                    continue;
                }
                let tid = match self.cur.reference_sequence_id() {
                    Some(t) => t,
                    None => continue,
                };
                let start = match self.cur.alignment_start() {
                    Some(p) => usize::from(p),
                    None => continue,
                };
                if self.current_tid != Some(tid) {
                    if self.current_tid.is_some() {
                        // Reference change: finish this one first.
                        self.held = true;
                        break;
                    }
                    self.current_tid = Some(tid);
                    self.current_pos = start;
                }
                if start > self.current_pos {
                    // Still ahead.
                    self.held = true;
                    break;
                }
                let buf = self.free.pop().unwrap_or_default();
                match decode_active(&self.cur, self.min_baseq, buf) {
                    Ok(ar) => self.active.push(ar),
                    Err(buf) => self.free.push(buf),
                }
            }

            self.retire_finished();

            col.clear();
            for a in &self.active {
                let offset = self.current_pos as i64 - a.ref_start as i64;
                if offset >= 0 && (offset as usize) < a.bases.len() {
                    let (b, q) = a.bases[offset as usize];
                    if b != b'*' {
                        col.push((b, q));
                    }
                }
            }

            let tid = match self.current_tid {
                Some(t) => t,
                None => return Ok(None),
            };
            let pos = self.current_pos;
            self.current_pos += 1;

            if !col.is_empty() {
                return Ok(Some((tid, pos)));
            }
            // No coverage: advance until something appears or EOF.
            if self.active.is_empty() && !self.held {
                if !self.fetch_next()? {
                    self.current_tid = None;
                    self.current_pos = 0;
                    continue;
                }
                self.held = true;
            }
        }
    }

    /// Drop reads that ended before the current position, keeping the survivors
    /// in order and reclaiming their buffers.
    fn retire_finished(&mut self) {
        let pos = self.current_pos;
        let mut keep = 0;
        for i in 0..self.active.len() {
            if self.active[i].end >= pos {
                self.active.swap(keep, i);
                keep += 1;
            }
        }
        for mut dead in self.active.drain(keep..) {
            dead.bases.clear();
            self.free.push(std::mem::take(&mut dead.bases));
        }
    }
}

/// Expand a record into per-reference-position `(base, qual)`, reusing `buf`.
/// Returns the buffer in `Err` when the record has no start.
fn decode_active(
    rec: &RecordBuf,
    min_baseq: u8,
    mut buf: Vec<(u8, u8)>,
) -> std::result::Result<ActiveRead, Vec<(u8, u8)>> {
    let Some(start) = rec.alignment_start() else {
        return Err(buf);
    };
    let start = usize::from(start);
    let seq = rec.sequence().as_ref();
    let quals = rec.quality_scores().as_ref();
    buf.clear();
    let mut read_idx = 0usize;
    let mut ref_idx = 0usize; // 0-based from `start`
    for op in rec.cigar().as_ref().iter() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    if read_idx >= seq.len() {
                        break;
                    }
                    let b = seq[read_idx];
                    let q = if read_idx < quals.len() {
                        quals[read_idx]
                    } else {
                        0
                    };
                    if q >= min_baseq {
                        buf.push((b, q));
                    } else {
                        buf.push((b'N', 0));
                    }
                    read_idx += 1;
                    ref_idx += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                read_idx += len;
            }
            Kind::Deletion | Kind::Skip => {
                for _ in 0..len {
                    buf.push((b'*', 0));
                    ref_idx += 1;
                }
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    let end = start + ref_idx.saturating_sub(1);
    Ok(ActiveRead {
        end,
        bases: buf,
        ref_start: start,
    })
}
