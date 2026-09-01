use std::fs::File;
use std::io::{self, BufReader, BufWriter, Cursor, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use memmap2::Mmap;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as AlignmentWrite;
use noodles_sam::header::record::value::{Map, map::Program, map::program::tag as pg_tag};

use crate::types::OutputFormat;

/// Resolve a memory hint into an absolute byte budget.
///
/// Accepts:
/// * `"auto"` / `""` — pick `max(min_bytes, frac_num/frac_den · total_RAM)`.
///   This is the default for the aligner-driven fused pipeline: on a 32 GB
///   box we end up with ~24 GB sort budget, which fits chr20 30× WGS in a
///   single chunk and skips the spill+external-merge phase entirely.
/// * `"768M"`, `"4G"`, `"512K"`, `"1234"` (bytes) — explicit size, matches
///   `samtools sort -m` conventions.
///
/// `frac_num`/`frac_den` express the auto fraction (e.g. 3/4 = 75 %). The
/// minimum floor (`min_bytes`) protects against `auto` picking absurdly
/// small numbers on memory-starved hosts.
pub fn resolve_memory_hint(
    spec: &str,
    min_bytes: usize,
    frac_num: u32,
    frac_den: u32,
) -> Result<usize> {
    let s = spec.trim();
    if s.is_empty() || s.eq_ignore_ascii_case("auto") {
        return Ok(auto_memory_bytes(min_bytes, frac_num, frac_den));
    }
    let (num_str, mult): (&str, usize) = if let Some(rest) = s.strip_suffix(['G', 'g']) {
        (rest, 1 << 30)
    } else if let Some(rest) = s.strip_suffix(['M', 'm']) {
        (rest, 1 << 20)
    } else if let Some(rest) = s.strip_suffix(['K', 'k']) {
        (rest, 1 << 10)
    } else {
        (s, 1)
    };
    let n: usize = num_str
        .trim()
        .parse()
        .with_context(|| format!("parse memory size `{spec}`"))?;
    Ok(n.saturating_mul(mult).max(1 << 20))
}

/// Total-RAM-relative buffer sizing.
///
/// Returns `max(min_bytes, total_ram * frac_num / frac_den)`.
///
/// We deliberately use *total* RAM rather than *available* RAM: the available
/// figure on Linux/Windows is conservative (subtracts caches that the kernel
/// would happily evict), and getting a 768M default when the box has 64 GB
/// of cache is exactly the failure mode we're trying to avoid.
pub fn auto_memory_bytes(min_bytes: usize, frac_num: u32, frac_den: u32) -> usize {
    use sysinfo::System;
    let mut sys = System::new();
    sys.refresh_memory();
    let total = sys.total_memory() as usize; // bytes since sysinfo 0.30
    if total == 0 {
        // sysinfo failed to detect; fall back to floor.
        return min_bytes;
    }
    let frac = total.saturating_mul(frac_num as usize) / (frac_den as usize).max(1);
    frac.max(min_bytes)
}

/// Program record info appended to header as a fresh `@PG` line.
/// Set `enable = false` to skip (mirrors samtools `--no-PG`).
#[derive(Clone, Debug)]
pub struct PgInfo {
    pub id_prefix: String,
    pub name: String,
    pub version: String,
    pub cli: String,
    pub enable: bool,
}

impl PgInfo {
    pub fn new(subcmd: &str, enable: bool) -> Self {
        // BAM header values forbid tabs/newlines; bwa-style escape so the original args stay readable.
        let cli = std::env::args()
            .map(|a| a.replace('\t', "\\t").replace('\n', "\\n"))
            .collect::<Vec<_>>()
            .join(" ");
        Self {
            id_prefix: format!("kira-bam.{subcmd}"),
            name: "kira-bam".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            cli,
            enable,
        }
    }
}

/// Append one `@PG` record to the header, chaining to the existing leaf via
/// `PP`.
///
/// Deliberately not `Programs::add`. noodles resolves an ID collision by
/// appending the *previous program's ID*, so re-running a recipe over its own
/// output grew IDs without bound
/// (`kira-bam.sort-kira-bam.sort-kira-bam.sort`); and when the header holds
/// several chains — which every merge produces — it adds one record per leaf,
/// duplicating this invocation across the header. samtools adds exactly one
/// record, with `ID`, `ID-1`, `ID-2`… on collision, so that is what this does.
pub fn append_pg(header: &mut Header, pg: &PgInfo) -> Result<()> {
    if !pg.enable {
        return Ok(());
    }
    let id = unique_pg_id(header, &pg.id_prefix);
    let previous = last_leaf_program_id(header);

    let mut builder = Map::<Program>::builder()
        .insert(pg_tag::NAME, pg.name.as_bytes())
        .insert(pg_tag::VERSION, pg.version.as_bytes())
        .insert(pg_tag::COMMAND_LINE, pg.cli.as_bytes());
    if let Some(pp) = &previous {
        builder = builder.insert(pg_tag::PREVIOUS_PROGRAM_ID, pp.as_slice());
    }
    let map = builder.build().context("build @PG record")?;

    header.programs_mut().as_mut().insert(id, map);
    Ok(())
}

/// First free ID in the `prefix`, `prefix-1`, `prefix-2`, … sequence.
fn unique_pg_id(header: &Header, prefix: &str) -> bstr::BString {
    use bstr::ByteSlice;

    let programs = header.programs().as_ref();
    let mut id = bstr::BString::from(prefix.as_bytes());
    let mut n = 1usize;
    while programs.contains_key(id.as_bstr()) {
        id = bstr::BString::from(format!("{prefix}-{n}").into_bytes());
        n += 1;
    }
    id
}

/// The last program that nothing else links back to — the end of the chain a
/// new record should extend. Headers can hold several chains (a merge unions
/// the inputs'), and the new record joins only the last of them, so one
/// invocation stays one `@PG` line.
fn last_leaf_program_id(header: &Header) -> Option<bstr::BString> {
    use bstr::ByteSlice;
    use std::collections::HashSet;

    let programs = header.programs().as_ref();
    let referenced: HashSet<&bstr::BStr> = programs
        .values()
        .filter_map(|m| m.other_fields().get(&pg_tag::PREVIOUS_PROGRAM_ID))
        .map(|v| v.as_bstr())
        .collect();
    programs
        .keys()
        .rfind(|id| !referenced.contains(id.as_bstr()))
        .cloned()
}

const BGZF_MAGIC: [u8; 2] = [0x1f, 0x8b];
const CRAM_MAGIC: [u8; 4] = *b"CRAM";

/// BGZF end-of-file marker: an empty deflate block, the last 28 bytes of a
/// well-formed BGZF file.
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// BAM's uncompressed magic, the first four bytes inside the BGZF stream.
const BAM_MAGIC: [u8; 4] = *b"BAM\x01";

/// Is this gzip member BGZF? BGZF is gzip with `FLG.FEXTRA` set and a `BC`
/// extra subfield carrying the block size; `gzip(1)` sets neither. Checking the
/// header structure separates the two without decompressing, exactly as
/// htslib's `hts_detect_format` does.
fn is_bgzf_header(head: &[u8]) -> bool {
    const FEXTRA: u8 = 0x04;
    if head.len() < 12 || head[3] & FEXTRA == 0 {
        return false;
    }
    let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
    let Some(extra) = head.get(12..12 + xlen) else {
        return false;
    };
    // Extra subfields are `SI1 SI2 SLEN(u16le) data[SLEN]`, packed back to back.
    let mut i = 0usize;
    while i + 4 <= extra.len() {
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        if extra[i] == b'B' && extra[i + 1] == b'C' {
            return true;
        }
        i += 4 + slen;
    }
    false
}

/// The first four decompressed bytes of a BGZF file, used to tell BAM from a
/// bgzipped SAM. A short or undecodable stream just reports "not BAM" — the
/// real parse error then comes from the reader, with its own context.
fn bgzf_leading_bytes(path: &Path) -> [u8; 4] {
    let mut out = [0u8; 4];
    if let Ok(file) = File::open(path) {
        let mut r = bgzf::io::Reader::new(file);
        // A partial read leaves the tail zeroed, which cannot match BAM_MAGIC.
        let _ = r.read(&mut out);
    }
    out
}

fn detect_format<P: AsRef<Path>>(path: P) -> io::Result<InputFormat> {
    let path = path.as_ref();
    let mut f = File::open(path)?;
    // Enough for the gzip header plus its extra field; BGZF's is 6 bytes.
    let mut buf = [0u8; 64];
    let n = f.read(&mut buf)?;
    let head = &buf[..n];

    let detected = if n >= 4 && head[..4] == CRAM_MAGIC {
        InputFormat::Cram
    } else if n >= 2 && head[..2] == BGZF_MAGIC {
        if !is_bgzf_header(head) {
            // Plain gzip. Feeding it to the BAM parser used to fail with an
            // opaque decode error several layers down; say what is actually
            // wrong and how to fix it.
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{} is plain gzip, not BGZF. kira-bam reads BGZF (what `bgzip` and every \
                     BAM writer produce); recompress with `bgzip`, or decompress first \
                     (`gunzip -c file | kira-bam ...`).",
                    path.display()
                ),
            ));
        }
        // BGZF carries BAM, but also bgzipped SAM — tell them apart by magic
        // rather than assuming, which sent every bgzipped SAM into the BAM
        // parser.
        if bgzf_leading_bytes(path) == BAM_MAGIC {
            InputFormat::Bam
        } else {
            InputFormat::BgzfSam
        }
    } else {
        InputFormat::Sam
    };

    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        let ext_l = ext.to_ascii_lowercase();
        let mismatch = matches!(
            (detected, ext_l.as_str()),
            (InputFormat::Bam, "sam")
                | (InputFormat::Bam, "cram")
                | (InputFormat::Sam, "bam")
                | (InputFormat::Sam, "cram")
                | (InputFormat::BgzfSam, "bam")
                | (InputFormat::BgzfSam, "cram")
                | (InputFormat::Cram, "sam")
                | (InputFormat::Cram, "bam")
        );
        if mismatch {
            eprintln!(
                "[kira-bam] warning: file {} has .{ext_l} extension but content looks like {:?}",
                path.display(),
                detected
            );
        }
    }
    Ok(detected)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InputFormat {
    Bam,
    Sam,
    /// SAM text inside a BGZF stream, as `bgzip file.sam` produces.
    BgzfSam,
    Cram,
}

/// Write-side IO options. Default = single-threaded BGZF compression, which
/// matches the legacy behaviour. Setting `compression_workers >= 2` switches
/// BAM writes to `bgzf::MultithreadedWriter`.
///
/// Recommended values:
/// * `0` or `1` — single-threaded (default).
/// * `min(threads - 1, 8)` — best general-case throughput for chr-scale
///   output, leaving 1 core for the producer thread and capping at 8 since
///   BGZF block compression saturates around there on consumer SSDs.
#[derive(Clone, Copy, Debug, Default)]
pub struct WriteOptions {
    /// Number of parallel BGZF compression workers. 0/1 = single-threaded.
    pub compression_workers: usize,
    /// BGZF compression level. `None` = noodles default (6, samtools-compatible).
    pub compression_level: Option<bgzf::io::writer::CompressionLevel>,
}

/// Read-side IO options. Default = stream the file via buffered `read()`.
/// Setting `mmap = true` switches BAM/SAM input to mmap-backed reads.
///
/// Mmap is a win when:
/// * the file is large enough that the kernel hasn't cached the whole thing
///   (one mmap + page faults vs many `read()` syscalls), or
/// * the file is hot in cache (the kernel turns reads into trivial page
///   table lookups; we still save the syscall overhead).
///
/// It's a wash or slightly worse on small files and on filesystems where
/// mmap requires more setup than the read path (e.g. some network mounts).
#[derive(Clone, Copy, Debug, Default)]
pub struct OpenOptions {
    /// Read via `memmap2::Mmap` instead of `File`. Has no effect on CRAM
    /// (CRAM's slice cache already does its own buffering). Has no effect
    /// on already-decompressed streams.
    pub mmap: bool,
}

pub struct BamReader {
    inner: ReaderImpl,
    header: Header,
    /// Reused by [`BamReader::read_summary`] on BAM input.
    raw: bam::Record,
    /// Reused by `read_summary` on SAM/CRAM.
    summary_fallback: RecordBuf,
}

/// The fixed-position record fields. See [`BamReader::read_summary`].
#[derive(Clone, Copy, Debug, Default)]
pub struct RecordSummary {
    pub flags: u16,
    /// `None` when MAPQ is the 255 "unavailable" sentinel.
    pub mapping_quality: Option<u8>,
    pub reference_sequence_id: Option<usize>,
    pub mate_reference_sequence_id: Option<usize>,
    pub alignment_start: Option<usize>,
    pub template_length: i32,
}

impl RecordSummary {
    fn fill_from_raw(&mut self, rec: &bam::Record) -> Result<()> {
        self.flags = u16::from(rec.flags());
        self.mapping_quality = rec.mapping_quality().map(u8::from);
        self.reference_sequence_id = rec.reference_sequence_id().transpose()?;
        self.mate_reference_sequence_id = rec.mate_reference_sequence_id().transpose()?;
        self.alignment_start = rec.alignment_start().transpose()?.map(usize::from);
        self.template_length = rec.template_length();
        Ok(())
    }

    fn fill_from_buf(&mut self, rec: &RecordBuf) {
        self.flags = u16::from(rec.flags());
        self.mapping_quality = rec.mapping_quality().map(u8::from);
        self.reference_sequence_id = rec.reference_sequence_id();
        self.mate_reference_sequence_id = rec.mate_reference_sequence_id();
        self.alignment_start = rec.alignment_start().map(usize::from);
        self.template_length = rec.template_length();
    }

    /// Reference id as samtools prints it: `-1` for unplaced.
    pub fn tid(&self) -> i32 {
        self.reference_sequence_id.map(|i| i as i32).unwrap_or(-1)
    }

    pub fn mate_tid(&self) -> i32 {
        self.mate_reference_sequence_id
            .map(|i| i as i32)
            .unwrap_or(-1)
    }
}

enum ReaderImpl {
    Bam(bam::io::Reader<bgzf::io::Reader<File>>),
    BamMmap(bam::io::Reader<bgzf::io::Reader<Cursor<MmapView>>>),
    Sam(sam::io::Reader<BufReader<File>>),
    /// SAM backed by an mmap'd file. The cursor consumes the byte slice
    /// directly; one syscall to mmap vs O(file_size / block_size) read()
    /// calls. Wins ~10-15 % on the aligner→sort handoff (tmp SAM is on
    /// disk, but the kernel already has it warm in page cache).
    SamMmap(sam::io::Reader<Cursor<MmapView>>),
    /// SAM text inside a BGZF stream (`bgzip file.sam`). mmap buys nothing here
    /// — the cost is decompression, not the read syscalls.
    SamBgzf(sam::io::Reader<BufReader<bgzf::io::Reader<File>>>),
    Cram(CramSource),
}

/// Container-at-a-time CRAM decoding. The smallest readable unit is a container
/// of thousands of records, and `Reader::records()` buffers one internally — so
/// the iterator must outlive each record. It borrows the reader, which a struct
/// cannot also hold, so this drives `read_container` directly.
struct CramSource {
    reader: cram::io::Reader<BufReader<File>>,
    container: cram::io::reader::Container,
    /// Empty by default, which only decodes CRAMs storing bases verbatim.
    repository: fasta::Repository,
    decoded: Vec<RecordBuf>,
    cursor: usize,
}

impl CramSource {
    fn read_record_buf(&mut self, header: &Header, out: &mut RecordBuf) -> Result<bool> {
        while self.cursor >= self.decoded.len() {
            if !self.fill(header)? {
                return Ok(false);
            }
        }
        *out = std::mem::take(&mut self.decoded[self.cursor]);
        self.cursor += 1;
        Ok(true)
    }

    /// Decode containers into `decoded` until one yields records. `false` at EOF.
    fn fill(&mut self, header: &Header) -> Result<bool> {
        self.decoded.clear();
        self.cursor = 0;
        loop {
            let n = self
                .reader
                .read_container(&mut self.container)
                .context("read CRAM container")?;
            if n == 0 {
                return Ok(false);
            }
            let compression_header = self
                .container
                .compression_header()
                .context("read CRAM compression header")?;
            for slice in self.container.slices() {
                let slice = slice.context("read CRAM slice")?;
                let (core_data_src, external_data_srcs) =
                    slice.decode_blocks().context("decode CRAM slice blocks")?;
                let records = slice
                    .records(
                        self.repository.clone(),
                        header,
                        &compression_header,
                        &core_data_src,
                        &external_data_srcs,
                    )
                    .context("decode CRAM slice records")?;
                for record in &records {
                    self.decoded.push(
                        RecordBuf::try_from_alignment_record(header, record)
                            .context("convert CRAM record")?,
                    );
                }
            }
            if !self.decoded.is_empty() {
                return Ok(true);
            }
        }
    }
}

/// `Mmap` newtype that exposes a `&[u8]` view via `Cursor`. Holds the mapping
/// alive for the reader's lifetime; dropping `MmapView` unmaps.
pub struct MmapView {
    _mmap: Mmap,
    ptr: *const u8,
    len: usize,
}

// SAFETY: the data is read-only, owned by `_mmap`, and we never expose `&mut`.
unsafe impl Send for MmapView {}
unsafe impl Sync for MmapView {}

impl AsRef<[u8]> for MmapView {
    fn as_ref(&self) -> &[u8] {
        // SAFETY: `ptr`/`len` come from `_mmap.as_ref()`, which is valid for
        // the lifetime of `_mmap`, and `_mmap` lives as long as `self`.
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
}

impl MmapView {
    fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open {} for mmap", path.display()))?;
        let mmap = unsafe { Mmap::map(&file).with_context(|| format!("mmap {}", path.display()))? };
        let bytes: &[u8] = mmap.as_ref();
        let ptr = bytes.as_ptr();
        let len = bytes.len();
        Ok(Self {
            _mmap: mmap,
            ptr,
            len,
        })
    }
}

impl BamReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::open_with_reference::<P, &Path>(path, None)
    }

    /// Open with optional reference FASTA (CRAM may require it for slice decode).
    pub fn open_with_reference<P, R>(path: P, reference: Option<R>) -> Result<Self>
    where
        P: AsRef<Path>,
        R: AsRef<Path>,
    {
        Self::open_with_options(path, reference, OpenOptions::default())
    }

    /// Open with full control over IO backing.
    ///
    /// Defaults match `open_with_reference`. Pass `OpenOptions { mmap: true, .. }`
    /// to read the file via mmap — wins on hot files (tmp SAM in the
    /// aligner→sort handoff is almost always in page cache by the time we
    /// open it, so mmap is one syscall vs N small `read()`s).
    pub fn open_with_options<P, R>(
        path: P,
        reference: Option<R>,
        options: OpenOptions,
    ) -> Result<Self>
    where
        P: AsRef<Path>,
        R: AsRef<Path>,
    {
        let path = path.as_ref();
        let fmt = detect_format(path).context("detect input format")?;
        match fmt {
            InputFormat::Bam => {
                if options.mmap {
                    let view = MmapView::open(path)?;
                    let cur = Cursor::new(view);
                    let bgz = bgzf::io::Reader::new(cur);
                    let mut reader = bam::io::Reader::from(bgz);
                    let header = reader.read_header().context("read BAM header")?;
                    return Ok(Self {
                        inner: ReaderImpl::BamMmap(reader),
                        header,
                        raw: bam::Record::default(),
                        summary_fallback: RecordBuf::default(),
                    });
                }
                let file = File::open(path).context("open BAM")?;
                let mut reader = bam::io::Reader::new(file);
                let header = reader.read_header().context("read BAM header")?;
                Ok(Self {
                    inner: ReaderImpl::Bam(reader),
                    header,
                    raw: bam::Record::default(),
                    summary_fallback: RecordBuf::default(),
                })
            }
            InputFormat::Sam => {
                if options.mmap {
                    let view = MmapView::open(path)?;
                    let cur = Cursor::new(view);
                    let mut reader = sam::io::Reader::new(cur);
                    let header = reader.read_header().context("read SAM header")?;
                    return Ok(Self {
                        inner: ReaderImpl::SamMmap(reader),
                        header,
                        raw: bam::Record::default(),
                        summary_fallback: RecordBuf::default(),
                    });
                }
                let file = File::open(path).context("open SAM")?;
                let mut reader = sam::io::Reader::new(BufReader::with_capacity(1 << 20, file));
                let header = reader.read_header().context("read SAM header")?;
                Ok(Self {
                    inner: ReaderImpl::Sam(reader),
                    header,
                    raw: bam::Record::default(),
                    summary_fallback: RecordBuf::default(),
                })
            }
            InputFormat::BgzfSam => {
                let file = File::open(path).context("open BGZF SAM")?;
                let bgz = bgzf::io::Reader::new(file);
                let mut reader = sam::io::Reader::new(BufReader::with_capacity(1 << 20, bgz));
                let header = reader.read_header().context("read SAM header (BGZF)")?;
                Ok(Self {
                    inner: ReaderImpl::SamBgzf(reader),
                    header,
                    raw: bam::Record::default(),
                    summary_fallback: RecordBuf::default(),
                })
            }
            InputFormat::Cram => {
                let file = File::open(path).context("open CRAM")?;
                let buf_file = BufReader::with_capacity(1 << 20, file);
                let mut builder = cram::io::reader::Builder::default();
                let repository = match reference {
                    Some(refpath) => build_fasta_repo(refpath.as_ref())?,
                    None => fasta::Repository::default(),
                };
                builder = builder.set_reference_sequence_repository(repository.clone());
                let mut reader = builder.build_from_reader(buf_file);
                // CRAM read_header() validates the magic, parses the file
                // definition, and returns the SAM header in one go.
                let header = reader.read_header().context("read CRAM header")?;
                Ok(Self {
                    inner: ReaderImpl::Cram(CramSource {
                        reader,
                        container: cram::io::reader::Container::default(),
                        repository,
                        decoded: Vec::new(),
                        cursor: 0,
                    }),
                    header,
                    raw: bam::Record::default(),
                    summary_fallback: RecordBuf::default(),
                })
            }
        }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn read_record_buf(&mut self, buf: &mut RecordBuf) -> Result<bool> {
        match &mut self.inner {
            ReaderImpl::Bam(r) => {
                let n = r
                    .read_record_buf(&self.header, buf)
                    .context("read BAM record")?;
                Ok(n > 0)
            }
            ReaderImpl::BamMmap(r) => {
                let n = r
                    .read_record_buf(&self.header, buf)
                    .context("read BAM record (mmap)")?;
                Ok(n > 0)
            }
            ReaderImpl::Sam(r) => {
                let n = r
                    .read_record_buf(&self.header, buf)
                    .context("read SAM record")?;
                Ok(n > 0)
            }
            ReaderImpl::SamMmap(r) => {
                let n = r
                    .read_record_buf(&self.header, buf)
                    .context("read SAM record (mmap)")?;
                Ok(n > 0)
            }
            ReaderImpl::SamBgzf(r) => {
                let n = r
                    .read_record_buf(&self.header, buf)
                    .context("read SAM record (BGZF)")?;
                Ok(n > 0)
            }
            ReaderImpl::Cram(r) => r.read_record_buf(&self.header, buf),
        }
    }

    /// Read only the fixed-position fields of the next record, skipping the
    /// name, CIGAR, sequence, qualities and tags that [`Self::read_record_buf`]
    /// would decode. Allocation-free for BAM; SAM and CRAM fall back to a full
    /// decode, having no lazy form.
    pub fn read_summary(&mut self, out: &mut RecordSummary) -> Result<bool> {
        match &mut self.inner {
            ReaderImpl::Bam(r) => {
                let n = r
                    .read_record(&mut self.raw)
                    .context("read BAM record header")?;
                if n == 0 {
                    return Ok(false);
                }
                out.fill_from_raw(&self.raw)?;
                Ok(true)
            }
            ReaderImpl::BamMmap(r) => {
                let n = r
                    .read_record(&mut self.raw)
                    .context("read BAM record header (mmap)")?;
                if n == 0 {
                    return Ok(false);
                }
                out.fill_from_raw(&self.raw)?;
                Ok(true)
            }
            _ => {
                let mut buf = std::mem::take(&mut self.summary_fallback);
                let more = self.read_record_buf(&mut buf)?;
                if more {
                    out.fill_from_buf(&buf);
                }
                self.summary_fallback = buf;
                Ok(more)
            }
        }
    }

    /// Byte offset into the *compressed* stream, for BAM input only.
    ///
    /// Lets a caller extrapolate a record count from a short sample — sizing a
    /// hash map, say — instead of paying a full decompression pass to count
    /// exactly.
    pub fn compressed_position(&self) -> Option<u64> {
        match &self.inner {
            ReaderImpl::Bam(r) => Some(r.get_ref().virtual_position().compressed()),
            ReaderImpl::BamMmap(r) => Some(r.get_ref().virtual_position().compressed()),
            _ => None,
        }
    }

    pub fn input_format(&self) -> InputFormat {
        match self.inner {
            ReaderImpl::Bam(_) | ReaderImpl::BamMmap(_) => InputFormat::Bam,
            ReaderImpl::Sam(_) | ReaderImpl::SamMmap(_) => InputFormat::Sam,
            ReaderImpl::SamBgzf(_) => InputFormat::BgzfSam,
            ReaderImpl::Cram(_) => InputFormat::Cram,
        }
    }
}

fn build_fasta_repo(path: &Path) -> Result<fasta::repository::Repository> {
    let adapter = fasta::repository::adapters::IndexedReader::new(
        fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("open FASTA {}", path.display()))?,
    );
    Ok(fasta::repository::Repository::new(adapter))
}

pub struct BamWriter {
    inner: WriterImpl,
    header: Header,
    header_written: bool,
}

enum WriterImpl {
    Bam(bam::io::Writer<bgzf::io::Writer<BufWriter<Box<dyn Write + Send>>>>),
    /// Multi-threaded BGZF compression — wraps the same `BufWriter` sink but
    /// fans out per-block deflate jobs over a worker pool. This is the
    /// default for the fused aligner→sort pipeline because BGZF compression
    /// dominates the write phase on chr-scale data (we observed 25-40 % of
    /// total sort wall in single-threaded compression).
    BamMt(bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<Box<dyn Write + Send>>>>),
    Sam(sam::io::Writer<BufWriter<Box<dyn Write + Send>>>),
    Cram(cram::io::Writer<BufWriter<Box<dyn Write + Send>>>),
}

impl BamWriter {
    pub fn create<P: AsRef<Path>>(
        path: Option<P>,
        header: Header,
        fmt: OutputFormat,
    ) -> Result<Self> {
        Self::create_with_reference::<P, &Path>(path, header, fmt, None)
    }

    pub fn create_with_reference<P, R>(
        path: Option<P>,
        header: Header,
        fmt: OutputFormat,
        reference: Option<R>,
    ) -> Result<Self>
    where
        P: AsRef<Path>,
        R: AsRef<Path>,
    {
        let sink: Box<dyn Write + Send> = match path {
            Some(p) => Box::new(File::create(p.as_ref()).context("create output")?),
            None => Box::new(io::stdout()),
        };
        Self::from_writer_with_reference(sink, header, fmt, reference)
    }

    /// Like `create_with_reference` but with full write-side option control.
    pub fn create_with_options<P, R>(
        path: Option<P>,
        header: Header,
        fmt: OutputFormat,
        reference: Option<R>,
        opts: WriteOptions,
    ) -> Result<Self>
    where
        P: AsRef<Path>,
        R: AsRef<Path>,
    {
        let sink: Box<dyn Write + Send> = match path {
            Some(p) => Box::new(File::create(p.as_ref()).context("create output")?),
            None => Box::new(io::stdout()),
        };
        Self::from_writer_with_options(sink, header, fmt, reference, opts)
    }

    pub fn from_writer(
        sink: Box<dyn Write + Send>,
        header: Header,
        fmt: OutputFormat,
    ) -> Result<Self> {
        Self::from_writer_with_reference::<&Path>(sink, header, fmt, None)
    }

    pub fn from_writer_with_reference<R>(
        sink: Box<dyn Write + Send>,
        header: Header,
        fmt: OutputFormat,
        reference: Option<R>,
    ) -> Result<Self>
    where
        R: AsRef<Path>,
    {
        Self::from_writer_with_options::<R>(sink, header, fmt, reference, WriteOptions::default())
    }

    /// Build a writer with full I/O option control. The relevant knobs:
    ///
    /// * `compression_workers` — number of BGZF compression workers for
    ///   BAM output. `0` selects single-threaded (current legacy default);
    ///   `>= 2` switches to `bgzf::MultithreadedWriter`. CRAM and SAM
    ///   outputs ignore this — CRAM has its own slice-level parallelism,
    ///   SAM is uncompressed.
    pub fn from_writer_with_options<R>(
        sink: Box<dyn Write + Send>,
        header: Header,
        fmt: OutputFormat,
        reference: Option<R>,
        opts: WriteOptions,
    ) -> Result<Self>
    where
        R: AsRef<Path>,
    {
        let buf = BufWriter::with_capacity(1 << 20, sink);
        let inner = match fmt {
            OutputFormat::Bam => {
                if opts.compression_workers >= 2 {
                    let workers = std::num::NonZero::new(opts.compression_workers)
                        .unwrap_or(std::num::NonZero::<usize>::MIN);
                    let mut builder = bgzf::io::multithreaded_writer::Builder::default()
                        .set_worker_count(workers);
                    if let Some(level) = opts.compression_level {
                        builder = builder.set_compression_level(level);
                    }
                    let bgz = builder.build_from_writer(buf);
                    WriterImpl::BamMt(bam::io::Writer::from(bgz))
                } else {
                    let mut builder = bgzf::io::writer::Builder::default();
                    if let Some(level) = opts.compression_level {
                        builder = builder.set_compression_level(level);
                    }
                    let bgz = builder.build_from_writer(buf);
                    WriterImpl::Bam(bam::io::Writer::from(bgz))
                }
            }
            OutputFormat::UncompressedBam => {
                // Actually uncompressed, which is what `view -u` and the sort
                // spill files ask for. This used to fall back to FAST to dodge
                // a noodles-bgzf level-0 problem; `uncompressed_bam_round_trips`
                // holds that fix in place.
                let bgz = bgzf::io::writer::Builder::default()
                    .set_compression_level(bgzf::io::writer::CompressionLevel::NONE)
                    .build_from_writer(buf);
                WriterImpl::Bam(bam::io::Writer::from(bgz))
            }
            OutputFormat::Sam => WriterImpl::Sam(sam::io::Writer::new(buf)),
            OutputFormat::Cram => {
                let mut builder = cram::io::writer::Builder::default();
                if let Some(refpath) = reference {
                    let repo = build_fasta_repo(refpath.as_ref())?;
                    builder = builder.set_reference_sequence_repository(repo);
                }
                let writer = builder.build_from_writer(buf);
                WriterImpl::Cram(writer)
            }
        };
        Ok(Self {
            inner,
            header,
            header_written: false,
        })
    }

    pub fn write_header(&mut self) -> Result<()> {
        if self.header_written {
            return Ok(());
        }
        match &mut self.inner {
            WriterImpl::Bam(w) => w.write_header(&self.header).context("write BAM header")?,
            WriterImpl::BamMt(w) => w.write_header(&self.header).context("write BAM header")?,
            WriterImpl::Sam(w) => w.write_header(&self.header).context("write SAM header")?,
            WriterImpl::Cram(w) => w.write_header(&self.header).context("write CRAM header")?,
        }
        self.header_written = true;
        Ok(())
    }

    pub fn write_record_buf(&mut self, record: &RecordBuf) -> Result<()> {
        // Binary formats are unparseable without a header — auto-emit if forgotten.
        if !self.header_written
            && matches!(
                self.inner,
                WriterImpl::Bam(_) | WriterImpl::BamMt(_) | WriterImpl::Cram(_)
            )
        {
            self.write_header()?;
        }
        match &mut self.inner {
            WriterImpl::Bam(w) => w
                .write_alignment_record(&self.header, record)
                .context("write BAM record")?,
            WriterImpl::BamMt(w) => w
                .write_alignment_record(&self.header, record)
                .context("write BAM record (mt)")?,
            WriterImpl::Sam(w) => w
                .write_alignment_record(&self.header, record)
                .context("write SAM record")?,
            WriterImpl::Cram(w) => w
                .write_alignment_record(&self.header, record)
                .context("write CRAM record")?,
        }
        Ok(())
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Stream pre-encoded BAM record bytes (from [`encode_records_into`]) through BGZF. BAM only.
    pub fn write_preencoded(&mut self, raw: &[u8]) -> Result<()> {
        if !self.header_written {
            self.write_header()?;
        }
        match &mut self.inner {
            WriterImpl::Bam(w) => w
                .get_mut()
                .write_all(raw)
                .context("write pre-encoded BAM")?,
            WriterImpl::BamMt(w) => w
                .get_mut()
                .write_all(raw)
                .context("write pre-encoded BAM (mt)")?,
            _ => anyhow::bail!("write_preencoded supports BAM output only"),
        }
        Ok(())
    }

    /// Finalise the output: EOF marker, then flush to the sink.
    ///
    /// `AlignmentWrite::finish` is a no-op for BAM and SAM in noodles; the EOF
    /// block and the flushes otherwise happen in `Drop`, which discards their
    /// errors and turns a full disk into a truncated file with exit status 0.
    pub fn finish(self) -> Result<()> {
        let Self { inner, header, .. } = self;
        match inner {
            WriterImpl::Bam(w) => {
                let mut sink = w
                    .into_inner()
                    .finish()
                    .context("write BGZF EOF block to BAM output")?;
                sink.flush().context("flush BAM output")?;
            }
            WriterImpl::BamMt(w) => {
                let mut mt = w.into_inner();
                let mut sink = mt.finish().context("write BGZF EOF block to BAM output")?;
                sink.flush().context("flush BAM output")?;
            }
            WriterImpl::Sam(w) => {
                w.into_inner().flush().context("flush SAM output")?;
            }
            WriterImpl::Cram(mut w) => {
                AlignmentWrite::finish(&mut w, &header).context("write CRAM EOF container")?;
                w.into_inner().flush().context("flush CRAM output")?;
            }
        }
        Ok(())
    }
}

/// Does this write failure concern the single record rather than the stream?
/// Callers skip invalid records but must abort on a failing sink, so anything
/// that is not explicitly record validity counts as fatal.
pub fn is_record_level_error(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|c| c.downcast_ref::<io::Error>())
        .any(|e| {
            matches!(
                e.kind(),
                io::ErrorKind::InvalidInput | io::ErrorKind::InvalidData
            )
        })
}

/// Is this a closed downstream pipe (`… | head`)? Callers exit successfully.
pub fn is_broken_pipe(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|c| c.downcast_ref::<io::Error>())
        .any(|e| e.kind() == io::ErrorKind::BrokenPipe)
}

/// Write one record, tolerating only the record-validity failures noodles
/// reports more strictly than samtools. Anything else — a full disk, a closed
/// pipe, a failing sink — propagates, so a truncated output can never exit 0.
///
/// `skipped` accumulates the tolerated failures; pass it to [`report_skipped`]
/// once the stream is done.
pub fn write_record_tolerant(
    writer: &mut BamWriter,
    record: &RecordBuf,
    tool: &str,
    skipped: &mut u64,
) -> Result<()> {
    if let Err(e) = writer.write_record_buf(record) {
        if !is_record_level_error(&e) {
            return Err(e);
        }
        *skipped += 1;
        if *skipped <= 5 {
            eprintln!("[kira-bam {tool}] skip record: {e:#}");
        } else if *skipped == 6 {
            eprintln!("[kira-bam {tool}] (further skips suppressed)");
        }
    }
    Ok(())
}

/// Summarise the records [`write_record_tolerant`] dropped.
pub fn report_skipped(tool: &str, skipped: u64) {
    if skipped > 0 {
        eprintln!("[kira-bam {tool}] {skipped} records skipped as invalid");
    }
}

/// Encode `records` into the raw BAM record stream consumed by
/// [`BamWriter::write_preencoded`], returning how many records noodles rejected.
///
/// Encoding targets a `Vec`, so a failure here is always record validity, never
/// I/O — the same tolerance the record-at-a-time paths apply.
pub fn encode_records_into(header: &Header, records: &[RecordBuf]) -> (Vec<u8>, u64) {
    let mut w = bam::io::Writer::from(Vec::<u8>::with_capacity(records.len() * 64));
    let mut skipped = 0u64;
    for r in records {
        if w.write_alignment_record(header, r).is_err() {
            skipped += 1;
        }
    }
    (w.into_inner(), skipped)
}

pub fn default_index_path(bam: &Path, csi: bool) -> PathBuf {
    let mut p = bam.to_path_buf().into_os_string();
    p.push(if csi { ".csi" } else { ".bai" });
    PathBuf::from(p)
}

pub fn install_thread_pool(threads: usize) {
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build_global();
}

pub use noodles_sam::alignment::RecordBuf as Record;

#[cfg(test)]
mod finish_tests {
    use super::*;

    /// Refuses every write, as a full disk would.
    struct FailingSink;

    impl Write for FailingSink {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("no space left on device"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Err(io::Error::other("no space left on device"))
        }
    }

    fn header() -> Header {
        Header::default()
    }

    /// Output is buffered, so nothing reaches the sink until finalisation —
    /// which is why discarding `finish`'s result truncated files silently.
    #[test]
    fn finish_reports_a_failing_sink() {
        for fmt in [
            OutputFormat::Bam,
            OutputFormat::UncompressedBam,
            OutputFormat::Sam,
        ] {
            let mut w = BamWriter::from_writer(Box::new(FailingSink), header(), fmt).unwrap();
            // Header writes land in the 1 MiB BufWriter, not on the sink.
            let _ = w.write_header();
            assert!(
                w.finish().is_err(),
                "finish() reported success for {fmt:?} despite a failing sink"
            );
        }
    }

    #[test]
    fn sink_failures_are_not_classified_as_record_level() {
        let stream = anyhow::Error::new(io::Error::other("no space left on device"))
            .context("write BAM record");
        assert!(!is_record_level_error(&stream));

        let broken = anyhow::Error::new(io::Error::new(io::ErrorKind::BrokenPipe, "pipe"))
            .context("write BAM record");
        assert!(!is_record_level_error(&broken));
        assert!(is_broken_pipe(&broken));

        let record =
            anyhow::Error::new(io::Error::new(io::ErrorKind::InvalidInput, "invalid CIGAR"))
                .context("write BAM record");
        assert!(is_record_level_error(&record));
        assert!(!is_broken_pipe(&record));
    }
}

#[cfg(test)]
mod uncompressed_tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::cigar::{Op, op::Kind};
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles_sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    /// `view -u` promises an uncompressed BAM, and `sort` writes its spill
    /// files that way to skip compression it never benefits from. Both were
    /// silently getting level 1 instead; this pins the level-0 round trip so
    /// the workaround cannot creep back unnoticed.
    #[test]
    fn uncompressed_bam_round_trips() {
        let len = NonZeroUsize::try_from(1000).unwrap();
        let header = Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(len))
            .build();

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("u.bam");

        let mut w =
            BamWriter::create(Some(&path), header.clone(), OutputFormat::UncompressedBam).unwrap();
        w.write_header().unwrap();
        for i in 0..100 {
            let rec = RecordBuf::builder()
                .set_name(format!("r{i}"))
                .set_flags(Flags::default())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(i + 1).unwrap())
                .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 4)]))
                .set_sequence(Sequence::from(b"ACGT".to_vec()))
                .set_quality_scores(QualityScores::from(vec![30u8; 4]))
                .build();
            w.write_record_buf(&rec).unwrap();
        }
        w.finish().unwrap();

        let mut reader = BamReader::open(&path).unwrap();
        assert_eq!(reader.header().reference_sequences().len(), 1);
        let mut rec = RecordBuf::default();
        let mut n = 0;
        while reader.read_record_buf(&mut rec).unwrap() {
            n += 1;
        }
        assert_eq!(n, 100, "records lost in the uncompressed round trip");
    }

    /// `bgzip file.sam` is a BGZF stream that is not a BAM. Assuming every
    /// gzip file was a BAM sent it to the BAM parser, which failed several
    /// layers down with nothing about the actual problem.
    #[test]
    fn bgzipped_sam_is_read_as_sam() {
        const SAM: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
             r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
             r2\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.sam.gz");
        let file = File::create(&path).unwrap();
        let mut w = bgzf::io::Writer::new(file);
        w.write_all(SAM.as_bytes()).unwrap();
        w.finish().unwrap();

        assert_eq!(detect_format(&path).unwrap(), InputFormat::BgzfSam);

        let mut reader = BamReader::open(&path).unwrap();
        assert_eq!(reader.header().reference_sequences().len(), 1);
        let mut rec = RecordBuf::default();
        let mut n = 0;
        while reader.read_record_buf(&mut rec).unwrap() {
            n += 1;
        }
        assert_eq!(n, 2);
    }

    /// Plain `gzip` output is not BGZF and cannot be read; the error has to say
    /// so rather than surface as a BAM decode failure.
    #[test]
    fn plain_gzip_is_rejected_by_name() {
        // A minimal gzip member: header with no FEXTRA, then an empty stored
        // deflate block. Only the header matters to detection.
        let mut bytes: Vec<u8> = vec![0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0x00, 0xff];
        bytes.extend_from_slice(&[0x03, 0x00]);
        bytes.extend_from_slice(&[0u8; 8]);

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("plain.sam.gz");
        std::fs::write(&path, &bytes).unwrap();

        let err = detect_format(&path).unwrap_err();
        assert!(
            err.to_string().contains("plain gzip"),
            "unexpected error: {err}"
        );
    }
}

#[cfg(test)]
mod pg_tests {
    use super::*;

    fn pg(subcmd: &str) -> PgInfo {
        PgInfo {
            id_prefix: format!("kira-bam.{subcmd}"),
            name: "kira-bam".to_string(),
            version: "0.0.0".to_string(),
            cli: "kira-bam test".to_string(),
            enable: true,
        }
    }

    fn ids(header: &Header) -> Vec<String> {
        header
            .programs()
            .as_ref()
            .keys()
            .map(|k| k.to_string())
            .collect()
    }

    fn pp_of(header: &Header, id: &str) -> Option<String> {
        header
            .programs()
            .as_ref()
            .get(id.as_bytes())
            .and_then(|m| m.other_fields().get(&pg_tag::PREVIOUS_PROGRAM_ID))
            .map(|v| v.to_string())
    }

    #[test]
    fn chains_through_pp() {
        let mut h = Header::default();
        append_pg(&mut h, &pg("sort")).unwrap();
        append_pg(&mut h, &pg("markdup")).unwrap();
        assert_eq!(ids(&h), ["kira-bam.sort", "kira-bam.markdup"]);
        assert_eq!(pp_of(&h, "kira-bam.sort"), None);
        assert_eq!(
            pp_of(&h, "kira-bam.markdup").as_deref(),
            Some("kira-bam.sort")
        );
    }

    /// Re-running the same subcommand over its own output must not collide, and
    /// must not grow the ID by the whole previous ID each time.
    #[test]
    fn repeated_runs_get_numbered_ids() {
        let mut h = Header::default();
        for _ in 0..4 {
            append_pg(&mut h, &pg("sort")).unwrap();
        }
        assert_eq!(
            ids(&h),
            [
                "kira-bam.sort",
                "kira-bam.sort-1",
                "kira-bam.sort-2",
                "kira-bam.sort-3"
            ]
        );
        assert_eq!(
            pp_of(&h, "kira-bam.sort-3").as_deref(),
            Some("kira-bam.sort-2")
        );
    }

    /// A merged header holds one chain per input. The new record joins the last
    /// one; it must not be duplicated across every chain.
    #[test]
    fn one_record_even_with_several_chains() {
        let mut h = Header::default();
        let leaf = |name: &str| {
            Map::<Program>::builder()
                .insert(pg_tag::NAME, name.as_bytes())
                .build()
                .unwrap()
        };
        h.programs_mut().as_mut().insert("bwa".into(), leaf("bwa"));
        h.programs_mut()
            .as_mut()
            .insert("bwa-1".into(), leaf("bwa-other"));

        append_pg(&mut h, &pg("merge")).unwrap();
        assert_eq!(ids(&h), ["bwa", "bwa-1", "kira-bam.merge"]);
        assert_eq!(pp_of(&h, "kira-bam.merge").as_deref(), Some("bwa-1"));
    }
}
