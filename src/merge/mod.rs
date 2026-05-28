use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

use crate::cli::MergeArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SortKind {
    Coordinate,
    QueryName,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct MergeKey {
    tid: u32,
    pos: u32,
    qname: Vec<u8>,
    flag: u16,
}

fn build_key(rec: &RecordBuf, _kind: SortKind) -> MergeKey {
    MergeKey {
        tid: rec
            .reference_sequence_id()
            .map(|i| i as u32)
            .unwrap_or(u32::MAX),
        pos: rec
            .alignment_start()
            .map(|p| usize::from(p) as u32)
            .unwrap_or(u32::MAX),
        qname: rec.name().map(|n| n.to_vec()).unwrap_or_default(),
        flag: u16::from(rec.flags()),
    }
}

impl Ord for MergeKey {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.tid, self.pos, self.flag, self.qname.as_slice()).cmp(&(
            other.tid,
            other.pos,
            other.flag,
            other.qname.as_slice(),
        ))
    }
}

impl PartialOrd for MergeKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct HeapItem {
    key: MergeKey,
    record: RecordBuf,
    src: usize,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.src == other.src
    }
}
impl Eq for HeapItem {}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .key
            .cmp(&self.key)
            .then_with(|| other.src.cmp(&self.src))
    }
}

pub fn run(args: MergeArgs) -> Result<()> {
    install_thread_pool(args.threads);

    if !args.force && args.output.exists() {
        anyhow::bail!(
            "output {} already exists; pass -f to overwrite",
            args.output.display()
        );
    }

    let kind = if args.name_sort {
        SortKind::QueryName
    } else {
        SortKind::Coordinate
    };

    let mut readers: Vec<BamReader> = args
        .inputs
        .iter()
        .map(|p: &PathBuf| {
            BamReader::open_with_reference(p, args.reference.as_deref())
                .with_context(|| format!("open {}", p.display()))
        })
        .collect::<Result<_>>()?;

    let mut merged_header = merge_headers(readers.iter().map(|r| r.header()).collect::<Vec<_>>())?;
    append_pg(&mut merged_header, &PgInfo::new("merge", !args.no_pg))?;

    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::with_capacity(readers.len());
    for (idx, r) in readers.iter_mut().enumerate() {
        let mut rec = RecordBuf::default();
        if r.read_record_buf(&mut rec)? {
            heap.push(HeapItem {
                key: build_key(&rec, kind),
                record: rec,
                src: idx,
            });
        }
    }

    let fmt = if args
        .output
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("cram"))
        .unwrap_or(false)
    {
        OutputFormat::Cram
    } else if args.uncompressed {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    };
    let mut writer = BamWriter::create_with_reference(
        Some(&args.output),
        merged_header,
        fmt,
        args.reference.as_deref(),
    )?;
    writer.write_header()?;

    while let Some(item) = heap.pop() {
        writer.write_record_buf(&item.record)?;
        let src = item.src;
        let mut rec = RecordBuf::default();
        if readers[src].read_record_buf(&mut rec)? {
            heap.push(HeapItem {
                key: build_key(&rec, kind),
                record: rec,
                src,
            });
        }
    }
    writer.finish()?;
    Ok(())
}

fn merge_headers(headers: Vec<&Header>) -> Result<Header> {
    use bstr::{BString, ByteSlice};

    let mut merged = headers
        .first()
        .copied()
        .cloned()
        .context("no input headers")?;
    for (i, h) in headers.iter().enumerate().skip(1) {
        // @SQ: must match by name + length. Reject mismatches loudly.
        if h.reference_sequences().len() != merged.reference_sequences().len() {
            anyhow::bail!(
                "input {i} has {} reference sequences, expected {}",
                h.reference_sequences().len(),
                merged.reference_sequences().len()
            );
        }
        for (name_h, sq_h) in h.reference_sequences() {
            match merged.reference_sequences().get(name_h.as_bstr()) {
                Some(sq_m) if sq_m.length() != sq_h.length() => {
                    anyhow::bail!(
                        "input {i}: @SQ {} has length {} vs expected {}",
                        name_h.as_bstr(),
                        usize::from(sq_h.length()),
                        usize::from(sq_m.length())
                    );
                }
                None => anyhow::bail!("input {i}: @SQ {} missing in base header", name_h.as_bstr()),
                _ => {}
            }
        }

        // @RG: dedup on ID, suffix duplicates with `-{idx}`.
        for (rg_id, rg_map) in h.read_groups() {
            let mut final_id: BString = rg_id.clone();
            let mut suffix = 1usize;
            while merged.read_groups().contains_key(final_id.as_bstr()) {
                final_id = BString::from(format!("{}-{suffix}", rg_id.as_bstr()));
                suffix += 1;
            }
            merged
                .read_groups_mut()
                .insert(final_id, rg_map.clone());
        }

        // @CO: append all comments.
        for c in h.comments() {
            merged.comments_mut().push(c.clone());
        }

        // @PG: append, but reuse noodles Programs::add so chains stay valid.
        for (pg_id, pg_map) in h.programs().as_ref() {
            let mut final_id: BString = pg_id.clone();
            let mut suffix = 1usize;
            while merged.programs().as_ref().contains_key(final_id.as_bstr()) {
                final_id = BString::from(format!("{}-{suffix}", pg_id.as_bstr()));
                suffix += 1;
            }
            // Use raw insert so we don't re-chain — preserve original PP links.
            merged
                .programs_mut()
                .as_mut()
                .insert(final_id, pg_map.clone());
        }
    }
    Ok(merged)
}
