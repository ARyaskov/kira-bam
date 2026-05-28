# kira-bam

High-performance BAM/SAM toolkit written in Rust 2024. Drop-in samtools-compatible
CLI (`view`/`sort`/`index`/`merge`/`markdup`/`flagstat`) plus a library API for
embedding directly into aligners and variant callers — designed to fuse the
classic `bwa | samtools sort | samtools markdup | samtools index` cascade into
one streaming pass.

## Goals

- 99% byte-compatible outputs against `samtools` on the same inputs.
- Wall-clock ≤ `samtools` on every command at the same `-@` thread count.
- Library-first design: `kira-ls-aligner` embeds `kira-bam` for fused
  align→sort→markdup→index, no intermediate SAM/BAM files.

## Installation

From source (Rust 1.95+):

```bash
cargo install --path .
```

Binary name: `kira-bam`.

## Commands

```bash
kira-bam view     [-bShH] [-f INT] [-F INT] [-q INT] [-c] [-@ N] [-o OUT] in.{sam,bam} [region ...]
kira-bam sort     [-n] [-u] [-@ N] [-m MEM] [-T DIR] [-o out.bam] in.bam
kira-bam index    [-b | -c] [--min-shift N] [--depth N] [-@ N] in.bam
kira-bam merge    [-n] [-u] [-f] [-@ N] -o out.bam in1.bam in2.bam ...
kira-bam markdup  [-r] [-d N] [-s STATS] [-u] [-@ N] in.sorted.bam out.bam
kira-bam flagstat [-@ N] in.bam
```

## Library use

```rust
use kira_bam::{BamReader, BamWriter};
use kira_bam::types::OutputFormat;

let mut reader = BamReader::open("in.bam")?;
let mut writer = BamWriter::create(Some("out.bam"), reader.header().clone(), OutputFormat::Bam)?;
writer.write_header()?;
let mut rec = kira_bam::io::Record::default();
while reader.read_record_buf(&mut rec)? {
    writer.write_record_buf(&rec)?;
}
writer.finish()?;
```

## Tests

`tests/golden.rs` does a per-file diff between `kira-bam` and `samtools` on a
directory of BAM fixtures. Skipped by default; enable with:

```bash
KIRA_BAM_SAMTOOLS=/path/to/samtools \
KIRA_BAM_FIXTURES=/path/to/bams \
cargo test --release
```

## License

MIT.
