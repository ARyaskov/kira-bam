//! Group records that share a read name, without holding the file in memory.
//!
//! Small inputs are grouped in a map and written straight out. Once the buffer
//! passes the memory budget, the run switches to samtools' approach: hash every
//! read name into one of `--buckets` temporary files, so all records of a name
//! land together, then group each bucket on its own. Peak memory becomes the
//! largest bucket rather than the whole BAM — the previous version kept every
//! record, cloned, and OOM'd on anything WGS-sized.

use std::collections::HashMap;
use std::hash::{BuildHasher, BuildHasherDefault};

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
use rustc_hash::FxHasher;

use crate::cli::CollateArgs;
use crate::io::{
    BamReader, BamWriter, PgInfo, append_pg, install_thread_pool, report_skipped,
    resolve_memory_hint, write_record_tolerant,
};
use crate::sort::approx_heap_size;
use crate::types::OutputFormat;

/// Floor for `-m auto`, matching the sort module's reasoning.
const AUTO_MIN_MEMORY_BYTES: usize = 512 * 1024 * 1024;
const AUTO_MEMORY_NUM: u32 = 3;
const AUTO_MEMORY_DEN: u32 = 4;

pub(crate) type Buckets = HashMap<Vec<u8>, Vec<RecordBuf>, BuildHasherDefault<FxHasher>>;

pub fn run(args: CollateArgs) -> Result<()> {
    install_thread_pool(args.threads);

    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let mut header = reader.header().clone();
    append_pg(&mut header, &PgInfo::new("collate", !args.no_pg))?;

    let budget = resolve_memory_hint(
        &args.memory,
        AUTO_MIN_MEMORY_BYTES,
        AUTO_MEMORY_NUM,
        AUTO_MEMORY_DEN,
    )?;
    let n_buckets = args.buckets.max(1);

    let mut out_path = args.output.clone();
    if out_path.extension().is_none() {
        out_path.set_extension("bam");
    }

    let mut buffered: Buckets = HashMap::default();
    let mut buffered_bytes = 0usize;
    let mut spill: Option<NameBuckets> = None;
    let mut rec = RecordBuf::default();
    let mut total = 0u64;

    while reader.read_record_buf(&mut rec)? {
        if args.fast {
            let f = u16::from(rec.flags());
            if (f & (0x100 | 0x800)) != 0 {
                continue;
            }
        }
        total += 1;

        if let Some(s) = spill.as_mut() {
            s.push(&rec)?;
            continue;
        }

        buffered_bytes += approx_heap_size(&rec);
        push_grouped(&mut buffered, std::mem::take(&mut rec));

        if buffered_bytes > budget {
            // Over budget: move everything buffered into bucket files and keep
            // streaming there.
            let mut s = NameBuckets::new(args.tmpdir.as_deref(), &header, n_buckets)?;
            for (_, recs) in buffered.drain() {
                for r in &recs {
                    s.push(r)?;
                }
            }
            buffered_bytes = 0;
            spill = Some(s);
        }
    }

    let mut writer = BamWriter::create(Some(&out_path), header.clone(), OutputFormat::Bam)?;
    writer.write_header()?;
    let mut skipped = 0u64;

    match spill {
        None => {
            for (_, recs) in buffered.drain() {
                for r in &recs {
                    write_record_tolerant(&mut writer, r, "collate", &mut skipped)?;
                }
            }
        }
        Some(s) => s.for_each_group(|recs| {
            for r in recs {
                write_record_tolerant(&mut writer, r, "collate", &mut skipped)?;
            }
            Ok(())
        })?,
    }

    writer.finish()?;
    report_skipped("collate", skipped);
    eprintln!(
        "[kira-bam collate] wrote {total} records to {}",
        out_path.display()
    );
    Ok(())
}

/// The temporary files records are hashed into, one per bucket.
///
/// Every record sharing a read name hashes to the same file, which is what
/// makes each bucket independently groupable in memory. Shared with `fastq`,
/// which needs the same "group by name without holding the file" guarantee.
pub(crate) struct NameBuckets {
    dir: tempfile::TempDir,
    writers: Vec<BamWriter>,
    paths: Vec<std::path::PathBuf>,
    hasher: BuildHasherDefault<FxHasher>,
}

impl NameBuckets {
    pub(crate) fn new(tmpdir: Option<&std::path::Path>, header: &Header, n: usize) -> Result<Self> {
        let dir = match tmpdir {
            Some(p) => tempfile::Builder::new()
                .prefix("kira-collate")
                .tempdir_in(p),
            None => tempfile::Builder::new().prefix("kira-collate").tempdir(),
        }
        .context("create collate temp dir")?;

        let mut writers = Vec::with_capacity(n);
        let mut paths = Vec::with_capacity(n);
        for i in 0..n {
            let path = dir.path().join(format!("bucket{i}.bam"));
            // Uncompressed: these files are read back once, on the same disk,
            // so compressing them only burns CPU.
            let mut w =
                BamWriter::create(Some(&path), header.clone(), OutputFormat::UncompressedBam)?;
            w.write_header()?;
            writers.push(w);
            paths.push(path);
        }
        Ok(Self {
            dir,
            writers,
            paths,
            hasher: BuildHasherDefault::default(),
        })
    }

    fn bucket_of(&self, name: &[u8]) -> usize {
        (self.hasher.hash_one(name) % self.writers.len() as u64) as usize
    }

    pub(crate) fn push(&mut self, rec: &RecordBuf) -> Result<()> {
        let name: &[u8] = rec.name().map(AsRef::as_ref).unwrap_or(&[]);
        let i = self.bucket_of(name);
        // A spill file that cannot be written is fatal: silently dropping
        // records here would quietly lose reads from the output.
        self.writers[i].write_record_buf(rec)
    }

    /// Read the buckets back one at a time, handing `f` each group of records
    /// that share a name. Only one bucket is resident at a time.
    pub(crate) fn for_each_group<F>(self, mut f: F) -> Result<()>
    where
        F: FnMut(&[RecordBuf]) -> Result<()>,
    {
        let Self {
            dir,
            writers,
            paths,
            ..
        } = self;
        for w in writers {
            w.finish()?;
        }
        for path in &paths {
            let mut reader = BamReader::open(path)
                .with_context(|| format!("reopen collate bucket {}", path.display()))?;
            let mut buckets: Buckets = HashMap::default();
            let mut rec = RecordBuf::default();
            while reader.read_record_buf(&mut rec)? {
                push_grouped(&mut buckets, std::mem::take(&mut rec));
            }
            for (_, recs) in buckets.drain() {
                f(&recs)?;
            }
        }
        // Held until every bucket has been read back; dropping it removes them.
        drop(dir);
        Ok(())
    }
}

/// File the record under its read name, allocating the key only for a name not
/// seen before.
pub(crate) fn push_grouped(buckets: &mut Buckets, rec: RecordBuf) {
    let name: &[u8] = rec.name().map(AsRef::as_ref).unwrap_or(&[]);
    match buckets.get_mut(name) {
        Some(v) => v.push(rec),
        None => {
            let key = name.to_vec();
            buckets.insert(key, vec![rec]);
        }
    }
}
