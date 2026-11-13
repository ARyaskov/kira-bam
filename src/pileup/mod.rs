//! Streaming pileup over a coordinate-sorted BAM.
//!
//! Yields `(tid, pos, [PileupObs])` in genomic order, keeping a sliding window
//! of active reads keyed by alignment end. A position's observations are one
//! per overlapping read, carrying the markup samtools mpileup renders:
//! strand (0x10), the read's MAPQ for the `^` marker at its first covered
//! position, an end flag for `$`, the base kind (aligned base, deletion `*`,
//! reference skip `>`/`<`) and an optional `+N<seq>` insertion and/or `-N`
//! deletion riding on the aligned base it follows, rendered insertion first.
//! Bases below `min_baseq` are not observed at
//! all — never an inflating `N` — but a position whose every observation was
//! filtered still yields an empty column, which mpileup renders as the
//! `0 * *` sentinel row. Reads with the samtools default exclusion
//! mask UNMAP|SECONDARY|QCFAIL|DUP (0x704; supplementary 0x800 retained) never
//! enter the window. A read's retained state is compact: read-level metadata
//! (MAPQ, strand, covered span) is stored once, aligned bases and qualities
//! are kept one byte each (bounded by the read, not the reference span), and
//! each `D`/`N` span is a single template (kind, quality, indel markup) plus a
//! length, expanded into a transient [`PileupObs`] only at the column that
//! covers it — a 50 kb reference skip costs one run, not 50 000 observations.
//! Nothing allocates per position in steady state: the column goes into a
//! caller-owned buffer, decode buffers are recycled as reads leave the window
//! (bounded by an explicit free-pool policy), and one record buffer is reused
//! throughout.
//!
//! SEQ is stored reference-forward in BAM, so no complementing is done: only
//! case and punctuation change with strand, exactly like samtools.

use std::sync::Arc;

use anyhow::Result;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::op::Kind;

use crate::io::RecordSource;

/// One observation in a pileup column: one read at one reference position,
/// with everything the samtools read-bases/qualities rendering needs.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct PileupObs {
    /// Base from SEQ, reference-forward. For deletions this is `b'*'`; for
    /// reference skips `b'>'`/`b'<'` by strand. `0` is the filtered sentinel:
    /// the slot keeps per-reference offsets aligned but is not an observation.
    pub base: u8,
    pub qual: u8,
    /// The read's MAPQ, printed after `^` at the read's first covered
    /// position.
    pub mapq: u8,
    /// Reverse strand (flag 0x10): bases print lowercase, skips as `<`.
    pub reverse: bool,
    /// First covered position of this alignment: renders `^`+MAPQ.
    pub is_start: bool,
    /// Last covered position of this alignment: renders `$`.
    pub is_end: bool,
    /// What the read does at this position.
    pub kind: ObsKind,
    /// `+N<seq>` attached to this observation: the read's inserted SEQ bases,
    /// `None` almost everywhere. Shared with the read's retained run, so
    /// materializing a column is a cheap `Arc` clone, not a per-column copy.
    pub ins: Option<Arc<[u8]>>,
    /// `-N` attached to this observation: `N` deleted reference bases follow
    /// the insertion text — a compound I+D on one entry renders
    /// `+N<ins>-N<del>`, insertion first, exactly like htslib's single
    /// per-entry indel field. `None` almost everywhere.
    pub del_len: Option<usize>,
}

/// What a read does at a reference position.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum ObsKind {
    /// An aligned base (M/=/X) or the filtered sentinel slot.
    #[default]
    Base,
    /// A deleted reference base: `*`.
    Deletion,
    /// A CIGAR N reference skip: `>` forward, `<` reverse.
    RefSkip,
}

pub type PileupColumn = (usize, usize, Vec<PileupObs>);

/// Bounds on the decode-buffer free pool (see [`recycle`]). A
/// deep/high-coverage region retires many reads at once; without a cap the
/// pool would pin every retired read's buffers for the rest of the run.
const FREE_MAX_COUNT: usize = 64;
const FREE_MAX_BUF_BYTES: usize = 1 << 20; // 1 MiB per pooled buffer

/// A compact reference-consuming run: one M/=/X, D or N span. Long spans are
/// stored as a single template plus length and expanded into a transient
/// [`PileupObs`] only at the column that covers them.
#[derive(Debug)]
struct Run {
    /// 0-based reference offset from the read's start where this run begins.
    start: usize,
    /// Reference bases this run consumes.
    len: usize,
    kind: RunKind,
    /// `+N`/`-N` indel markup riding on this run's LAST observation (the
    /// base the following I/D op attaches to). `None` almost everywhere.
    end_ins: Option<Arc<[u8]>>,
    end_del_len: Option<usize>,
}

#[derive(Debug)]
enum RunKind {
    /// M/=/X: aligned bases at `seq`/`qual` offsets `base_off .. base_off+len`.
    Match { base_off: usize },
    /// D: `*` entries carrying the quality of the base after the deletion.
    Del { qual: u8 },
    /// N: `>`/`<` entries carrying the quality of the base after the skip.
    Skip { qual: u8 },
}

#[derive(Debug, Default)]
pub struct ActiveRead {
    pub end: usize,
    /// 1-based reference start.
    pub ref_start: usize,
    /// Read-level metadata, stored once per read instead of once per position.
    mapq: u8,
    reverse: bool,
    /// Total reference-consuming length (M/=/X/D/N) — the read's covered
    /// span — for the `$` marker and retirement.
    ref_len: usize,
    /// Full SEQ in query order, including inserted and clipped bases.
    seq: Vec<u8>,
    /// Quality scores in query order; empty when the record has none.
    qual: Vec<u8>,
    /// Compact reference-consuming runs in order.
    runs: Vec<Run>,
}

impl ActiveRead {
    /// The transient observation at reference offset `off` (0-based from the
    /// read's start), or `None` when no run covers it or the base fails `-Q`.
    fn observe(&self, off: usize, min_baseq: u8) -> Option<PileupObs> {
        let i = self.runs.partition_point(|r| r.start <= off);
        let run = self.runs.get(i.checked_sub(1)?)?;
        if off >= run.start + run.len {
            // A SEQ shorter than the CIGAR's M span leaves no run here.
            return None;
        }
        let last = off + 1 == run.start + run.len;
        let (base, qual, kind) = match run.kind {
            RunKind::Match { base_off } => {
                let qi = base_off + (off - run.start);
                let b = *self.seq.get(qi)?;
                let q = self.qual.get(qi).copied().unwrap_or(0);
                (if q < min_baseq { 0 } else { b }, q, ObsKind::Base)
            }
            RunKind::Del { qual } => (
                if qual < min_baseq { 0 } else { b'*' },
                qual,
                ObsKind::Deletion,
            ),
            RunKind::Skip { qual } => {
                let b = if self.reverse { b'<' } else { b'>' };
                (if qual < min_baseq { 0 } else { b }, qual, ObsKind::RefSkip)
            }
        };
        if base == 0 {
            return None; // filtered: not an observation
        }
        Some(PileupObs {
            base,
            qual,
            mapq: self.mapq,
            reverse: self.reverse,
            is_start: off == 0,
            is_end: off + 1 == self.ref_len,
            kind,
            ins: if last { run.end_ins.clone() } else { None },
            del_len: if last { run.end_del_len } else { None },
        })
    }
}

pub struct PileupIter {
    reader: Box<dyn RecordSource>,
    active: Vec<ActiveRead>,
    /// Decode buffers of reads that left the window.
    free: Vec<ActiveRead>,
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
    pub fn new(reader: impl RecordSource + 'static, min_mapq: u8, min_baseq: u8) -> Self {
        Self {
            reader: Box::new(reader),
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
    pub fn next_column_into(&mut self, col: &mut Vec<PileupObs>) -> Result<Option<(usize, usize)>> {
        loop {
            // Pull in records that start at or before current_pos for current_tid.
            // Exits either at EOF or having parked one record via `hold`.
            loop {
                if !self.fetch_next()? {
                    break;
                }
                let flags = u16::from(self.cur.flags());
                // samtools mpileup's default exclusion: UNMAP is always
                // dropped and `--ff` defaults to SECONDARY|QCFAIL|DUP, i.e.
                // the 0x704 mask; SUPPLEMENTARY (0x800) is retained. Hardcoded
                // in the shared engine — the flag set mpileup and consensus
                // are built around.
                if (flags & 0x704) != 0 {
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
                let mut ar = self.free.pop().unwrap_or_default();
                match decode_active(&self.cur, &mut ar) {
                    Ok(()) => self.active.push(ar),
                    Err(()) => self.free.push(ar),
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
                if offset >= 0 {
                    // Deletions (`*`) and reference skips (`>`/`<`) are
                    // observations too: they count toward depth and render in
                    // the column, exactly like htslib's pileup. Only the
                    // filtered sentinel is dropped.
                    if let Some(obs) = a.observe(offset as usize, self.min_baseq) {
                        col.push(obs);
                    }
                }
            }

            // Reads cover this position. When every observation failed `-Q`
            // the column stays empty and the caller renders samtools' `0 * *`
            // sentinel row — a genuinely uncovered position never reaches
            // here, because the window-drain jump above skips it. The window
            // shrinks as `pos` advances, so this cannot spin.
            return Ok(Some((tid, pos)));
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
        for dead in self.active.drain(keep..) {
            recycle(&mut self.free, dead);
        }
    }
}

/// Recycle a retired read's buffers for the next record, dropping the whole
/// buffer when the pool is full or the buffer is too large to pin.
/// shortcut: the byte cap uses each buffer's current *capacity* (its peak
/// read length), not the live length; revisit if reads longer than
/// `FREE_MAX_BUF_BYTES` become the common case.
fn recycle(free: &mut Vec<ActiveRead>, mut buf: ActiveRead) {
    let bytes =
        buf.seq.capacity() + buf.qual.capacity() + buf.runs.capacity() * std::mem::size_of::<Run>();
    if bytes > FREE_MAX_BUF_BYTES || free.len() >= FREE_MAX_COUNT {
        return; // drop: the Vecs (and any Arc ins bytes) are deallocated
    }
    buf.seq.clear();
    buf.qual.clear();
    buf.runs.clear();
    free.push(buf);
}

/// Expand a record into compact per-read state in `ar`, reusing its buffers.
/// Returns `Err(())` when the record has no alignment start (nothing filled).
fn decode_active(rec: &RecordBuf, ar: &mut ActiveRead) -> std::result::Result<(), ()> {
    let Some(start) = rec.alignment_start() else {
        return Err(());
    };
    let start = usize::from(start);
    let seq = rec.sequence().as_ref();
    let quals = rec.quality_scores().as_ref();
    let reverse = (u16::from(rec.flags()) & 0x10) != 0;
    // noodles maps the MAPQ-255 "unavailable" marker to `None`; samtools
    // renders it as the 126 cap (`~`), never as 0.
    let mapq = rec.mapping_quality().map(u8::from).unwrap_or(255);

    // Total reference-consuming length (M/=/X/D/N): the read's last covered
    // position, where the `$` marker goes.
    let mut ref_len = 0usize;
    for op in rec.cigar().as_ref().iter() {
        match op.kind() {
            Kind::Match
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
            | Kind::Deletion
            | Kind::Skip => ref_len += op.len(),
            _ => {}
        }
    }

    ar.seq.clear();
    ar.seq.extend_from_slice(seq);
    ar.qual.clear();
    ar.qual.extend_from_slice(quals);
    ar.runs.clear();

    let ops = rec.cigar().as_ref();
    let mut read_idx = 0usize;
    let mut ref_idx = 0usize; // 0-based from `start`
    let mut i = 0;
    // shortcut: CIGARs that start with D/I, and reads whose first or last
    // covered base fails `-Q`, keep the marker/indel placement implemented
    // here; it matches htslib's entry semantics but is unverified against a
    // differential oracle for those exact edges.
    while i < ops.len() {
        let kind = ops[i].kind();
        let mut len = ops[i].len();
        i += 1;
        // htslib normalizes consecutive same-type I/D (1I1I→2I, 1D2D→3D);
        // without the merge the second op would overwrite the first's markup.
        if kind == Kind::Insertion || kind == Kind::Deletion {
            while i < ops.len() && ops[i].kind() == kind {
                len += ops[i].len();
                i += 1;
            }
        }
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // A SEQ shorter than the CIGAR's M span truncates the run; the
                // tail of the span is left uncovered, exactly like the
                // per-base loop it replaces.
                let actual = if read_idx < seq.len() {
                    (seq.len() - read_idx).min(len)
                } else {
                    0
                };
                if actual > 0 {
                    ar.runs.push(Run {
                        start: ref_idx,
                        len: actual,
                        kind: RunKind::Match { base_off: read_idx },
                        end_ins: None,
                        end_del_len: None,
                    });
                }
                read_idx += actual;
                ref_idx += actual;
            }
            Kind::Insertion => {
                // `+N<seq>` rides on the last emitted observation; a missing
                // inserted base renders `N` like the per-base loop.
                if let Some(last) = ar.runs.last_mut() {
                    last.end_ins = Some(ins_arc(seq, read_idx, len));
                }
                read_idx += len;
            }
            Kind::SoftClip => {
                read_idx += len;
            }
            Kind::Deletion => {
                // `-N` rides on the observation immediately before the D op;
                // the deleted positions themselves are `*` observations
                // carrying the quality of the base AFTER the D.
                if let Some(last) = ar.runs.last_mut() {
                    last.end_del_len = Some(len);
                }
                let q = quals.get(read_idx).copied().unwrap_or(0);
                ar.runs.push(Run {
                    start: ref_idx,
                    len,
                    kind: RunKind::Del { qual: q },
                    end_ins: None,
                    end_del_len: None,
                });
                ref_idx += len;
            }
            Kind::Skip => {
                let q = quals.get(read_idx).copied().unwrap_or(0);
                ar.runs.push(Run {
                    start: ref_idx,
                    len,
                    kind: RunKind::Skip { qual: q },
                    end_ins: None,
                    end_del_len: None,
                });
                ref_idx += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    ar.mapq = mapq;
    ar.reverse = reverse;
    ar.ref_start = start;
    ar.ref_len = ref_len;
    ar.end = start + ref_len.saturating_sub(1);
    Ok(())
}

/// The inserted SEQ bases at `read_idx..read_idx+len`, `N` filler for bases
/// past the end of SEQ — shared as an `Arc` so column materialization clones
/// the handle, not the bytes.
fn ins_arc(seq: &[u8], read_idx: usize, len: usize) -> Arc<[u8]> {
    let mut v = Vec::with_capacity(len);
    for k in 0..len {
        v.push(seq.get(read_idx + k).copied().unwrap_or(b'N'));
    }
    Arc::from(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{BamReader, BamWriter, RecordSource};
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

    /// Collect every column the iterator yields to exhaustion on a worker
    /// thread; also available as raw columns.
    fn observed(records: &[RecordBuf], min_baseq: u8) -> Vec<PileupColumn> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("in.bam");
        write_bam(&path, records);

        let (tx, rx) = std::sync::mpsc::channel();
        let handle = std::thread::spawn(move || {
            let reader = BamReader::open(&path).unwrap();
            let mut it = PileupIter::new(reader, 0, min_baseq);
            let mut out = Vec::new();
            while let Some(col) = it.next_column().unwrap() {
                out.push(col);
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

    /// Drive the iterator to exhaustion on a worker thread and record only
    /// column positions and depth. A regression that stops it terminating
    /// would otherwise hang the whole test run instead of failing; here it
    /// fails in 30 s with a clear message.
    fn columns(records: &[RecordBuf]) -> Vec<(usize, usize, usize)> {
        observed(records, 0)
            .into_iter()
            .map(|(tid, pos, col)| (tid, pos, col.len()))
            .collect()
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

    /// A read with explicit per-base qualities.
    fn readq(name: &str, tid: usize, start: usize, quals: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(tid)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, quals.len())]))
            .set_sequence(Sequence::from(vec![b'A'; quals.len()]))
            .set_quality_scores(QualityScores::from(quals.to_vec()))
            .build()
    }

    /// A read with a deletion in the middle: `2M 1D <len-3>M`.
    fn read_del(name: &str, tid: usize, start: usize, quals: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(tid)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Match, quals.len() - 2),
            ]))
            .set_sequence(Sequence::from(vec![b'A'; quals.len()]))
            .set_quality_scores(QualityScores::from(quals.to_vec()))
            .build()
    }

    /// A `1M 1D 1M` read with explicit per-base qualities.
    fn read_del_mdm(name: &str, quals: &[u8]) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(10).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Match, 1),
            ]))
            .set_sequence(Sequence::from(b"AG".to_vec()))
            .set_quality_scores(QualityScores::from(quals.to_vec()))
            .build()
    }

    #[test]
    fn deletion_entry_carries_the_following_base_quality_and_passes_q() {
        // resolve_cigar2 keeps the query position on the base AFTER the D,
        // and bam_plcmd.c applies `-Q` to the `*` entry like every other:
        // quals [40, 10] with -Q 13 drop the deletion at position 11, while
        // the aligned bases (40) survive at 10 and 12.
        let cols = observed(&[read_del_mdm("r1", &[40, 10])], 13);
        let del = cols
            .iter()
            .find(|&&(_, p, _)| p == 11)
            .expect("the read still spans the deleted position");
        assert!(del.2.is_empty(), "qual 10 < 13 must drop the `*`: {cols:?}");
        // Same read without -Q: the `*` stays and carries the following
        // base's quality 10 (`K` = 10+33 in the quality column).
        let cols = observed(&[read_del_mdm("r1", &[40, 10])], 0);
        let del = cols.iter().find(|&&(_, p, _)| p == 11).unwrap();
        assert_eq!(del.2[0].base, b'*', "{cols:?}");
        assert_eq!(del.2[0].qual, 10, "following base's quality: {cols:?}");
    }

    #[test]
    fn trailing_deletion_has_no_following_base_so_zero_quality() {
        // A D at the very end of the read: htslib's `qpos < l_qseq` bound
        // yields quality 0, which `-Q 13` filters out; without -Q the `*`
        // entry survives with quality 0 (`!` in the quality column).
        let r = RecordBuf::builder()
            .set_name("r")
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(10).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Deletion, 1),
            ]))
            .set_sequence(Sequence::from(b"A".to_vec()))
            .set_quality_scores(QualityScores::from(vec![40]))
            .build();
        let cols = observed(std::slice::from_ref(&r), 13);
        let del = cols
            .iter()
            .find(|&&(_, p, _)| p == 11)
            .expect("the read still spans the deleted position");
        assert!(del.2.is_empty(), "qual 0 < 13 must drop the `*`: {cols:?}");
        let cols = observed(&[r], 0);
        let del = cols.iter().find(|&&(_, p, _)| p == 11).unwrap();
        assert_eq!(del.2[0].base, b'*', "{cols:?}");
        assert_eq!(
            del.2[0].qual, 0,
            "no following base means quality 0: {cols:?}"
        );
    }

    #[test]
    fn deletion_counts_as_depth_and_stays_in_the_column() {
        let cols = observed(&[read_del("r1", 0, 10, &[30; 5])], 0);
        // The deleted reference position (12) still yields a column and the
        // deletion is one observation, not a dropped read. It carries the
        // quality of the base following the D, like an htslib pileup entry
        // (resolve_cigar2 keeps the query position there).
        let del = cols
            .iter()
            .find(|&&(_, p, _)| p == 12)
            .expect("deletion column missing");
        assert_eq!(
            del.2,
            vec![PileupObs {
                base: b'*',
                qual: 30,
                mapq: 255,
                kind: ObsKind::Deletion,
                ..Default::default()
            }],
            "{cols:?}"
        );
        assert_eq!(cols.len(), 6, "one observation per ref position: {cols:?}");
    }

    /// A read jumping a reference skip: `2M 3N 3M` covering 10..=17.
    fn read_skip(name: &str, tid: usize, start: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::default())
            .set_reference_sequence_id(tid)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(Kind::Match, 2),
                Op::new(Kind::Skip, 3),
                Op::new(Kind::Match, 2),
            ]))
            .set_sequence(Sequence::from(b"AAAA".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 4]))
            .build()
    }

    #[test]
    fn base_below_min_baseq_is_not_an_n_or_a_depth_observation() {
        // Qual 5 with -Q 13: pos 10 keeps its read but every observation is
        // filtered, so it yields an EMPTY column (depth 0) — the caller
        // renders samtools' `0 * *` row — never `[(N, 1)]`. The filtered
        // first slot also takes the `^` marker with it, exactly like a
        // samtools entry failing `-Q`.
        let mut quals = vec![30u8; 5];
        quals[0] = 5;
        let cols = observed(&[readq("r1", 0, 10, &quals)], 13);
        assert_eq!(
            cols,
            vec![
                (0, 10, vec![]),
                (0, 11, vec![base_obs(b'A', 30)]),
                (0, 12, vec![base_obs(b'A', 30)]),
                (0, 13, vec![base_obs(b'A', 30)]),
                (
                    0,
                    14,
                    vec![PileupObs {
                        base: b'A',
                        qual: 30,
                        mapq: 255,
                        is_end: true,
                        ..Default::default()
                    }]
                ),
            ],
            "{cols:?}",
        );
    }

    #[test]
    fn reference_skip_is_an_observation_never_a_deletion() {
        // 2M covers 10-11, the 3N skip spans 12-14, 2M covers 15-16. The skip
        // is an observation like a deletion: it yields columns with `>`
        // (forward strand), counts toward depth, and carries the quality of
        // the base following the skip. `^` sits on the first covered position
        // (10) and `$` on the last (16).
        let cols = observed(&[read_skip("r1", 0, 10)], 0);
        assert_eq!(
            cols,
            vec![
                (
                    0,
                    10,
                    vec![PileupObs {
                        base: b'A',
                        qual: 30,
                        mapq: 255,
                        is_start: true,
                        ..Default::default()
                    }]
                ),
                (0, 11, vec![base_obs(b'A', 30)]),
                (
                    0,
                    12,
                    vec![PileupObs {
                        base: b'>',
                        qual: 30,
                        mapq: 255,
                        kind: ObsKind::RefSkip,
                        ..Default::default()
                    }]
                ),
                (
                    0,
                    13,
                    vec![PileupObs {
                        base: b'>',
                        qual: 30,
                        mapq: 255,
                        kind: ObsKind::RefSkip,
                        ..Default::default()
                    }]
                ),
                (
                    0,
                    14,
                    vec![PileupObs {
                        base: b'>',
                        qual: 30,
                        mapq: 255,
                        kind: ObsKind::RefSkip,
                        ..Default::default()
                    }]
                ),
                (0, 15, vec![base_obs(b'A', 30)]),
                (
                    0,
                    16,
                    vec![PileupObs {
                        base: b'A',
                        qual: 30,
                        mapq: 255,
                        is_end: true,
                        ..Default::default()
                    }]
                ),
            ],
            "the skip span must render `>` observations, one per position: {cols:?}",
        );
        assert!(
            !cols
                .iter()
                .any(|(_, _, c)| c.iter().any(|o| o.base == b'*')),
            "no skipped position may print as a deletion: {cols:?}",
        );
    }

    fn base_obs(base: u8, qual: u8) -> PileupObs {
        PileupObs {
            base,
            qual,
            // Reads without a MAPQ field carry the 255 "unavailable" marker.
            mapq: 255,
            ..Default::default()
        }
    }

    /// Reads carrying the samtools default exclusion mask must never reach
    /// the window; supplementary (0x800) must.
    #[test]
    fn default_excluded_flags_are_dropped_and_supplementary_is_kept() {
        use noodles_sam::alignment::record::Flags;
        let mk = |name: &str, flags: Flags| {
            RecordBuf::builder()
                .set_name(name)
                .set_flags(flags)
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(10).unwrap())
                .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5)]))
                .set_sequence(Sequence::from(b"ACGTA".to_vec()))
                .set_quality_scores(QualityScores::from(vec![30u8; 5]))
                .build()
        };
        let cols = columns(&[
            mk("sec", Flags::SECONDARY),
            mk("qcfail", Flags::QC_FAIL),
            mk("dup", Flags::DUPLICATE),
            mk("unmap", Flags::UNMAPPED),
            mk("supp", Flags::SUPPLEMENTARY),
        ]);
        assert_eq!(
            cols.len(),
            5,
            "only the supplementary read piles up: {cols:?}"
        );
        assert!(cols.iter().all(|&(_, _, d)| d == 1), "{cols:?}");
    }

    /// The decode-buffer free pool must not pin a huge read's buffers: a
    /// buffer above the byte cap is dropped, a normal one is pooled for reuse.
    #[test]
    fn free_pool_drops_an_oversized_buffer() {
        struct Empty;
        impl RecordSource for Empty {
            fn read_record_buf(&mut self, _out: &mut RecordBuf) -> Result<bool> {
                Ok(false)
            }
        }
        let mut it = PileupIter::new(Empty, 0, 0);

        let big = ActiveRead {
            seq: Vec::with_capacity(FREE_MAX_BUF_BYTES + 1),
            ..Default::default()
        };
        recycle(&mut it.free, big);
        assert!(
            it.free.is_empty(),
            "a buffer above the byte cap must be dropped, not pooled"
        );

        let small = ActiveRead {
            seq: vec![b'A'; 4],
            ..Default::default()
        };
        recycle(&mut it.free, small);
        assert_eq!(it.free.len(), 1, "a normal buffer is pooled for reuse");
    }
}
