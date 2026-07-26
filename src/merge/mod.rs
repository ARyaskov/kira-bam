use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::path::PathBuf;

use anyhow::{Context, Result};
use bstr::{BStr, BString};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

use crate::cli::MergeArgs;
use crate::io::{BamReader, BamWriter, PgInfo, append_pg, install_thread_pool};
// Reuse the sort module's samtools-compatible ordering so `merge` and `sort`
// agree on tie-breaks (coordinate: tid/pos/reverse-strand; name: strnum qname +
// READ1/READ2 bits), with equal-key records resolved by input/source order.
use crate::sort::{SortKey, SortKind, build_key};
use crate::types::OutputFormat;

struct HeapItem {
    key: SortKey,
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

    // Clone input headers up front: we need them to both build the merged header
    // and compute per-input tid remaps, but the readers are borrowed mutably for
    // the streaming merge below.
    let input_headers: Vec<Header> = readers.iter().map(|r| r.header().clone()).collect();
    let (mut merged_header, rg_renames) =
        merge_headers(input_headers.iter().collect::<Vec<_>>())?;
    append_pg(&mut merged_header, &PgInfo::new("merge", !args.no_pg))?;
    // Per-input map: local reference id → merged reference id (matched by @SQ name).
    let tid_maps = build_tid_maps(&input_headers, &merged_header);

    // A streaming k-way coordinate merge is only valid if every input is already
    // ordered consistently with the merged @SQ list — i.e. its local→merged tid
    // map is monotonic. If an input's @SQ order is a genuine permutation of the
    // merged order, remapping alone would emit records out of coordinate order, so
    // we refuse loudly rather than write a mis-sorted BAM. (Name-sorted merges key
    // on qname, so tid order is irrelevant there.)
    if kind == SortKind::Coordinate {
        for (i, map) in tid_maps.iter().enumerate() {
            if !is_sorted_compatible(map) {
                anyhow::bail!(
                    "input {} ({}) has @SQ in an order incompatible with the merged header; \
                     a coordinate merge needs all inputs ordered consistently. Re-sort it \
                     against the target header (`kira-bam sort`) or reheader it first.",
                    i,
                    args.inputs[i].display()
                );
            }
        }
    }

    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::with_capacity(readers.len());
    for (idx, r) in readers.iter_mut().enumerate() {
        let mut rec = RecordBuf::default();
        if r.read_record_buf(&mut rec)? {
            remap_ids(&mut rec, &tid_maps[idx]);
            remap_read_group(&mut rec, &rg_renames[idx]);
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
            let mut rec = item.record;
        if readers[src].read_record_buf(&mut rec)? {
            remap_ids(&mut rec, &tid_maps[src]);
            remap_read_group(&mut rec, &rg_renames[src]);
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

/// Per-input map from a read-group ID as stored in that input's records to the
/// one it was given in the merged header.
type RgRenames = HashMap<BString, BString>;

/// Build the merged header, plus the read-group renames it forced.
fn merge_headers(headers: Vec<&Header>) -> Result<(Header, Vec<RgRenames>)> {
    use bstr::ByteSlice;

    let mut merged = headers
        .first()
        .copied()
        .cloned()
        .context("no input headers")?;
    let mut rg_renames: Vec<RgRenames> = vec![RgRenames::default(); headers.len()];
    for (i, h) in headers.iter().enumerate().skip(1) {
        // @SQ: build the UNION across inputs, preserving first-seen order. A name
        // shared by several inputs must agree on length (hard error); names unique
        // to a later input are appended. Records from every input are then remapped
        // out of their local tid space into this merged space (see `build_tid_maps`
        // / `remap_ids`), so inputs whose @SQ lists are reordered — or are proper
        // subsets/supersets of each other — no longer silently mis-reference reads.
        for (name_h, sq_h) in h.reference_sequences() {
            // Take the length out by value first so we don't hold an immutable
            // borrow of `merged` across the `insert` in the `None` arm.
            let existing_len = merged
                .reference_sequences()
                .get(name_h.as_bstr())
                .map(|sq_m| sq_m.length());
            match existing_len {
                Some(len_m) => {
                    if len_m != sq_h.length() {
                        anyhow::bail!(
                            "input {i}: @SQ {} has length {} vs expected {}",
                            name_h.as_bstr(),
                            usize::from(sq_h.length()),
                            usize::from(len_m)
                        );
                    }
                }
                None => {
                    merged
                        .reference_sequences_mut()
                        .insert(name_h.clone(), sq_h.clone());
                }
            }
        }

        // An ID shared by two inputs is the same read group when the records
        // agree; otherwise the later one is suffixed and the rename recorded, or
        // its records would rebind to the first input's read group.
        for (rg_id, rg_map) in h.read_groups() {
            if merged.read_groups().get(rg_id.as_bstr()) == Some(rg_map) {
                continue;
            }
            let mut final_id: BString = rg_id.clone();
            let mut suffix = 1usize;
            while merged.read_groups().contains_key(final_id.as_bstr()) {
                final_id = BString::from(format!("{}-{suffix}", rg_id.as_bstr()));
                suffix += 1;
            }
            if final_id != *rg_id {
                rg_renames[i].insert(rg_id.clone(), final_id.clone());
            }
            merged.read_groups_mut().insert(final_id, rg_map.clone());
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
    Ok((merged, rg_renames))
}

/// For each input header, map its local reference id (the integer tid stored in
/// every record) to the corresponding id in the merged header, matched by @SQ
/// name. `None` means the name is absent from the merged header — impossible
/// when the merged header is the union of all inputs, but handled defensively.
fn build_tid_maps(inputs: &[Header], merged: &Header) -> Vec<Vec<Option<usize>>> {
    use std::collections::HashMap;
    let merged_idx: HashMap<String, usize> = merged
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect();
    inputs
        .iter()
        .map(|h| {
            h.reference_sequences()
                .iter()
                .map(|(name, _)| merged_idx.get(&name.to_string()).copied())
                .collect()
        })
        .collect()
}

/// True if remapping preserves coordinate order: the merged ids reached by
/// visiting local tids 0,1,2,… (the order a coordinate-sorted input emits them)
/// are non-decreasing. A `false` result means the input's @SQ order is a
/// permutation of the merged order and cannot be stream-merged.
fn is_sorted_compatible(map: &[Option<usize>]) -> bool {
    let mut last: Option<usize> = None;
    for &m in map {
        if let Some(v) = m {
            if let Some(l) = last
                && v < l
            {
                return false;
            }
            last = Some(v);
        }
    }
    true
}

/// Rewrite a record's reference id and mate reference id from one input's local
/// tid space into the merged tid space. Unmapped ends (`None` id) are left
/// alone; an id missing from the map is left unchanged (defensive — see
/// [`build_tid_maps`]).
fn remap_ids(rec: &mut RecordBuf, map: &[Option<usize>]) {
    if let Some(id) = rec.reference_sequence_id()
        && let Some(&Some(new)) = map.get(id)
    {
        *rec.reference_sequence_id_mut() = Some(new);
    }
    if let Some(mid) = rec.mate_reference_sequence_id()
        && let Some(&Some(new)) = map.get(mid)
    {
        *rec.mate_reference_sequence_id_mut() = Some(new);
    }
}

/// Point the record's `RG` tag at this input's read group under its merged name.
fn remap_read_group(rec: &mut RecordBuf, renames: &RgRenames) {
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;

    if renames.is_empty() {
        return;
    }
    let new = match rec.data().get(&Tag::READ_GROUP) {
        Some(Value::String(id)) => match renames.get(id.as_ref() as &BStr) {
            Some(new) => new.clone(),
            None => return,
        },
        _ => return,
    };
    rec.data_mut().insert(Tag::READ_GROUP, Value::String(new));
}
