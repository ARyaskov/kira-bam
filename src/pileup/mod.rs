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
    /// Where the held record starts. Lets an empty window jump straight to the
    /// next event — the next contig, or the far side of a coverage gap —
    /// instead of stepping one reference base at a time.
    held_tid: usize,
    held_pos: usize,
    current_tid: Option<usize>,
    current_pos: usize,
    /// Last record position seen, to verify the input is coordinate-sorted.
    last: Option<(usize, usize)>,
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
            held_tid: 0,
            held_pos: 0,
            current_tid: None,
            current_pos: 0,
            last: None,
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

    /// Park `cur` until the sweep reaches it, remembering where that is.
    fn hold(&mut self, tid: usize, start: usize) {
        self.held = true;
        self.held_tid = tid;
        self.held_pos = start;
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
            // Exits either at EOF or having parked one record via `hold`.
            loop {
                if !self.fetch_next()? {
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
                // The window only ever moves forward, so a record behind it is
                // dropped rather than piled up. Say so instead of quietly
                // reporting depths that are missing reads.
                if let Some((last_tid, last_start)) = self.last
                    && (tid < last_tid || (tid == last_tid && start < last_start))
                {
                    anyhow::bail!(
                        "input is not coordinate-sorted: reference {tid} position {start} \
                         follows reference {last_tid} position {last_start}. A pileup sweeps \
                         positions in order and would silently skip the out-of-order reads — \
                         sort it first (`kira-bam sort`)."
                    );
                }
                self.last = Some((tid, start));

                if self.current_tid.is_none() {
                    // First placed record decides where the sweep starts.
                    self.current_tid = Some(tid);
                    self.current_pos = start;
                }
                if self.current_tid != Some(tid) || start > self.current_pos {
                    // A later contig, or further along this one. Park it: the
                    // window still has to be drained where we are.
                    self.hold(tid, start);
                    break;
                }
                let buf = self.free.pop().unwrap_or_default();
                match decode_active(&self.cur, self.min_baseq, buf) {
                    Ok(ar) => self.active.push(ar),
                    Err(buf) => self.free.push(buf),
                }
            }

            self.retire_finished();

            if self.active.is_empty() {
                // Nothing covers this position and the window is drained, so no
                // column between here and the parked record can be non-empty:
                // jump to it. This is what closes a coverage gap in one step
                // instead of one reference base per iteration — and, when the
                // parked record is on another contig, what actually moves the
                // sweep onto it. Leaving `current_tid` behind here spun over
                // the rest of the old reference forever.
                if self.held {
                    self.current_tid = Some(self.held_tid);
                    self.current_pos = self.held_pos;
                    continue;
                }
                // Not holding anything means the fill loop stopped at EOF.
                return Ok(None);
            }

            let tid = match self.current_tid {
                Some(t) => t,
                None => return Ok(None),
            };
            let pos = self.current_pos;
            self.current_pos += 1;

            col.clear();
            for a in &self.active {
                let offset = pos as i64 - a.ref_start as i64;
                if offset >= 0 && (offset as usize) < a.bases.len() {
                    let (b, q) = a.bases[offset as usize];
                    if b != b'*' {
                        col.push((b, q));
                    }
                }
            }

            if !col.is_empty() {
                return Ok(Some((tid, pos)));
            }
            // Covered only by deletions/skips. The window shrinks as `pos`
            // advances, so this cannot spin.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{BamReader, BamWriter};
    use crate::types::OutputFormat;
    use noodles_core::Position;
    use noodles_sam::Header;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles_sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;
    use std::path::Path;

    fn header() -> Header {
        let len = NonZeroUsize::try_from(1_000_000).unwrap();
        Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(len))
            .add_reference_sequence("chr2", Map::<ReferenceSequence>::new(len))
            .build()
    }

    /// A 5 bp perfect match at `(tid, start)`.
    fn read(name: &str, tid: usize, start: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(tid)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5)]))
            .set_sequence(Sequence::from(b"ACGTA".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 5]))
            .build()
    }

    fn write_bam(path: &Path, records: &[RecordBuf]) {
        let mut w = BamWriter::create(Some(path), header(), OutputFormat::Bam).unwrap();
        w.write_header().unwrap();
        for r in records {
            w.write_record_buf(r).unwrap();
        }
        w.finish().unwrap();
    }

    /// Drive the iterator to exhaustion on a worker thread. A regression that
    /// stops it terminating would otherwise hang the whole test run instead of
    /// failing; here it fails in 30 s with a clear message.
    fn columns(records: &[RecordBuf]) -> Vec<(usize, usize, usize)> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("in.bam");
        write_bam(&path, records);

        let (tx, rx) = std::sync::mpsc::channel();
        let handle = std::thread::spawn(move || {
            let reader = BamReader::open(&path).unwrap();
            let mut it = PileupIter::new(reader, 0, 0);
            let mut out = Vec::new();
            while let Some((tid, pos, col)) = it.next_column().unwrap() {
                out.push((tid, pos, col.len()));
                if out.len() > 10_000 {
                    break;
                }
            }
            let _ = tx.send(out);
            // Keep `dir` alive until the reader is dropped.
            drop(dir);
        });
        let out = rx
            .recv_timeout(std::time::Duration::from_secs(30))
            .expect("PileupIter did not terminate");
        handle.join().unwrap();
        out
    }

    #[test]
    fn advances_to_the_next_contig() {
        // Two reads on different references. The sweep used to set `held` on
        // the contig change without ever adopting the new tid, so it spun on
        // chr1 forever and emitted nothing at all.
        let cols = columns(&[read("r1", 0, 10), read("r2", 1, 10)]);
        let tids: Vec<usize> = cols.iter().map(|c| c.0).collect();
        assert!(tids.contains(&0), "chr1 columns missing: {cols:?}");
        assert!(tids.contains(&1), "chr2 columns missing: {cols:?}");
        assert_eq!(cols.len(), 10, "expected 5 columns per read: {cols:?}");
    }

    #[test]
    fn covers_every_base_of_a_single_read() {
        let cols = columns(&[read("r1", 0, 10)]);
        assert_eq!(
            cols,
            vec![(0, 10, 1), (0, 11, 1), (0, 12, 1), (0, 13, 1), (0, 14, 1)]
        );
    }

    #[test]
    fn jumps_across_a_coverage_gap() {
        // Half a megabase of empty reference between the reads. Stepping one
        // position per iteration made this O(gap); it must be O(1).
        let cols = columns(&[read("r1", 0, 10), read("r2", 0, 500_000)]);
        assert_eq!(cols.len(), 10);
        assert_eq!(cols[4], (0, 14, 1));
        assert_eq!(cols[5], (0, 500_000, 1));
    }

    #[test]
    fn overlapping_reads_stack_in_one_column() {
        let cols = columns(&[read("r1", 0, 10), read("r2", 0, 12)]);
        // 10,11 = 1 read; 12,13,14 = 2 reads; 15,16 = 1 read.
        assert_eq!(
            cols,
            vec![
                (0, 10, 1),
                (0, 11, 1),
                (0, 12, 2),
                (0, 13, 2),
                (0, 14, 2),
                (0, 15, 1),
                (0, 16, 1)
            ]
        );
    }

    #[test]
    fn empty_input_terminates() {
        assert!(columns(&[]).is_empty());
    }

    #[test]
    fn unsorted_input_is_an_error_not_silently_dropped_reads() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("unsorted.bam");
        // The second read is behind the first, so the sweep has already passed
        // it by the time it arrives.
        write_bam(&path, &[read("r1", 0, 100), read("r2", 0, 10)]);

        let reader = BamReader::open(&path).unwrap();
        let mut it = PileupIter::new(reader, 0, 0);
        let mut err = None;
        for _ in 0..100 {
            match it.next_column() {
                Ok(Some(_)) => {}
                Ok(None) => break,
                Err(e) => {
                    err = Some(e);
                    break;
                }
            }
        }
        let err = err.expect("out-of-order input should be reported");
        assert!(
            err.to_string().contains("not coordinate-sorted"),
            "unexpected error: {err:#}"
        );
    }
}
