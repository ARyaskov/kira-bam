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
pub fn resolve_memory_hint(spec: &str, min_bytes: usize, frac_num: u32, frac_den: u32) -> Result<usize> {
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

/// Append a `@PG` record to the header, chaining to the existing leaf via `PP`.
pub fn append_pg(header: &mut Header, pg: &PgInfo) -> Result<()> {
    if !pg.enable {
        return Ok(());
    }
    let map = Map::<Program>::builder()
        .insert(pg_tag::NAME, pg.name.as_bytes())
        .insert(pg_tag::VERSION, pg.version.as_bytes())
        .insert(pg_tag::COMMAND_LINE, pg.cli.as_bytes())
        .build()
        .context("build @PG record")?;
    header
        .programs_mut()
        .add(pg.id_prefix.as_bytes(), map)
        .context("attach @PG to header chain")?;
    Ok(())
}

const BGZF_MAGIC: [u8; 2] = [0x1f, 0x8b];
const CRAM_MAGIC: [u8; 4] = *b"CRAM";

fn detect_format<P: AsRef<Path>>(path: P) -> io::Result<InputFormat> {
    let path = path.as_ref();
    let mut f = File::open(path)?;
    let mut buf = [0u8; 4];
    let n = f.read(&mut buf)?;
    let detected = if n >= 4 && buf == CRAM_MAGIC {
        InputFormat::Cram
    } else if n >= 2 && buf[..2] == BGZF_MAGIC {
        InputFormat::Bam
    } else {
        InputFormat::Sam
    };
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        let ext_l = ext.to_ascii_lowercase();
        let mismatch = match (detected, ext_l.as_str()) {
            (InputFormat::Bam, "sam") | (InputFormat::Bam, "cram") => true,
            (InputFormat::Sam, "bam") | (InputFormat::Sam, "cram") => true,
            (InputFormat::Cram, "sam") | (InputFormat::Cram, "bam") => true,
            _ => false,
        };
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
    Cram(cram::io::Reader<BufReader<File>>),
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
        let mmap = unsafe {
            Mmap::map(&file).with_context(|| format!("mmap {}", path.display()))?
        };
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
                    });
                }
                let file = File::open(path).context("open BAM")?;
                let mut reader = bam::io::Reader::new(file);
                let header = reader.read_header().context("read BAM header")?;
                Ok(Self {
                    inner: ReaderImpl::Bam(reader),
                    header,
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
                    });
                }
                let file = File::open(path).context("open SAM")?;
                let mut reader = sam::io::Reader::new(BufReader::with_capacity(1 << 20, file));
                let header = reader.read_header().context("read SAM header")?;
                Ok(Self {
                    inner: ReaderImpl::Sam(reader),
                    header,
                })
            }
            InputFormat::Cram => {
                let file = File::open(path).context("open CRAM")?;
                let buf_file = BufReader::with_capacity(1 << 20, file);
                let mut builder = cram::io::reader::Builder::default();
                if let Some(refpath) = reference {
                    let repo = build_fasta_repo(refpath.as_ref())?;
                    builder = builder.set_reference_sequence_repository(repo);
                }
                let mut reader = builder.build_from_reader(buf_file);
                // CRAM read_header() validates the magic, parses the file
                // definition, and returns the SAM header in one go.
                let header = reader.read_header().context("read CRAM header")?;
                Ok(Self {
                    inner: ReaderImpl::Cram(reader),
                    header,
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
            ReaderImpl::Cram(r) => {
                // CRAM doesn't have an in-place read_record_buf; use the iterator and try_from.
                let mut iter = r.records(&self.header);
                match iter.next() {
                    None => Ok(false),
                    Some(result) => {
                        let rec = result.context("read CRAM record")?;
                        *buf = RecordBuf::try_from_alignment_record(&self.header, &rec)
                            .context("convert CRAM record")?;
                        Ok(true)
                    }
                }
            }
        }
    }

    pub fn input_format(&self) -> InputFormat {
        match self.inner {
            ReaderImpl::Bam(_) | ReaderImpl::BamMmap(_) => InputFormat::Bam,
            ReaderImpl::Sam(_) | ReaderImpl::SamMmap(_) => InputFormat::Sam,
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
                    let bgz = bgzf::io::MultithreadedWriter::with_worker_count(workers, buf);
                    WriterImpl::BamMt(bam::io::Writer::from(bgz))
                } else {
                    let bgz = bgzf::io::Writer::new(buf);
                    WriterImpl::Bam(bam::io::Writer::from(bgz))
                }
            }
            OutputFormat::UncompressedBam => {
                // FAST not NONE — works around noodles-bgzf 0.47 level-0 round-trip regression.
                let bgz = bgzf::io::writer::Builder::default()
                    .set_compression_level(bgzf::io::writer::CompressionLevel::FAST)
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

    pub fn finish(mut self) -> Result<()> {
        match &mut self.inner {
            WriterImpl::Bam(w) => {
                let _ = AlignmentWrite::finish(w, &self.header);
            }
            WriterImpl::BamMt(w) => {
                let _ = AlignmentWrite::finish(w, &self.header);
            }
            WriterImpl::Sam(w) => {
                let _ = AlignmentWrite::finish(w, &self.header);
            }
            WriterImpl::Cram(w) => {
                let _ = AlignmentWrite::finish(w, &self.header);
            }
        }
        Ok(())
    }
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
