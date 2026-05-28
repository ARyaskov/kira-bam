//! Streaming pileup over a coordinate-sorted BAM.
//!
//! Yields, in genomic order, `(tid, pos, depth, [(base, qual)])`. Maintains a
//! sliding window of active reads keyed by their alignment end so we never look
//! back at records already past their range.

use anyhow::Result;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::RecordBuf;

use crate::io::BamReader;

#[derive(Debug)]
pub struct ActiveRead {
    pub end: usize,
    /// Pre-decoded per-reference-position (base, qual) for the window we have not yet emitted.
    pub bases: Vec<(u8, u8)>,
    /// Index into `bases` for the next ref position to consume.
    pub cursor: usize,
    /// 1-based reference start.
    pub ref_start: usize,
}

pub struct PileupIter {
    reader: BamReader,
    active: Vec<ActiveRead>,
    pending: Option<RecordBuf>,
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
            pending: None,
            current_tid: None,
            current_pos: 0,
            min_mapq,
            min_baseq,
        }
    }

    /// Reads the next record from the BAM, or pops a pending one. Returns true if we
    /// pulled in a record (so the caller knows there's more work).
    fn fetch_next(&mut self) -> Result<Option<RecordBuf>> {
        if let Some(r) = self.pending.take() {
            return Ok(Some(r));
        }
        let mut rec = RecordBuf::default();
        if self.reader.read_record_buf(&mut rec)? {
            Ok(Some(rec))
        } else {
            Ok(None)
        }
    }

    /// Return the next pileup column. Returns Ok(None) at EOF.
    pub fn next_column(&mut self) -> Result<Option<(usize, usize, Vec<(u8, u8)>)>> {
        loop {
            // Pull in records that start at or before current_pos for current_tid.
            loop {
                let rec = self.fetch_next()?;
                let Some(rec) = rec else {
                    if self.active.is_empty() && self.current_tid.is_none() {
                        return Ok(None);
                    }
                    break;
                };
                let flags = u16::from(rec.flags());
                if (flags & 0x4) != 0 {
                    continue;
                }
                if let Some(mq) = rec.mapping_quality() {
                    if u8::from(mq) < self.min_mapq {
                        continue;
                    }
                }
                let tid = match rec.reference_sequence_id() {
                    Some(t) => t,
                    None => continue,
                };
                let start = match rec.alignment_start() {
                    Some(p) => usize::from(p),
                    None => continue,
                };
                if self.current_tid != Some(tid) {
                    if self.current_tid.is_some() {
                        // Reference change. Hold record, finish current ref first.
                        self.pending = Some(rec);
                        break;
                    }
                    self.current_tid = Some(tid);
                    self.current_pos = start;
                }
                if start > self.current_pos {
                    // Still ahead; save for later.
                    self.pending = Some(rec);
                    break;
                }
                // Eligible: decode bases and add to active.
                if let Some(ar) = decode_active(&rec, self.min_baseq) {
                    self.active.push(ar);
                }
            }

            // Drop reads that have ended.
            self.active.retain(|a| a.end >= self.current_pos);

            // Build column for current position.
            let mut col: Vec<(u8, u8)> = Vec::new();
            for a in &mut self.active {
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
                return Ok(Some((tid, pos, col)));
            }
            // No coverage: advance until we find something or EOF.
            if self.active.is_empty() && self.pending.is_none() {
                // Try to load next batch.
                if self.fetch_next()?.is_none() {
                    self.current_tid = None;
                    self.current_pos = 0;
                    continue;
                }
            }
        }
    }
}

fn decode_active(rec: &RecordBuf, min_baseq: u8) -> Option<ActiveRead> {
    let start = usize::from(rec.alignment_start()?);
    let seq = rec.sequence().as_ref();
    let quals = rec.quality_scores().as_ref();
    let mut bases: Vec<(u8, u8)> = Vec::new();
    let mut read_idx = 0usize;
    let mut ref_idx = 0usize; // 0-based from `start`
    for op in rec.cigar().as_ref().iter() {
        let len = op.len() as usize;
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
                        bases.push((b, q));
                    } else {
                        bases.push((b'N', 0));
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
                    bases.push((b'*', 0));
                    ref_idx += 1;
                }
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    let end = start + ref_idx.saturating_sub(1);
    Some(ActiveRead {
        end,
        bases,
        cursor: 0,
        ref_start: start,
    })
}
