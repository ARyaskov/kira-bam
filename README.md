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

Unit tests (including a `proptest` that pins the AVX2 quality-sum to its scalar
reference across the full byte domain) run with a plain `cargo test`.

`tests/golden.rs` cross-checks `kira-bam` against a real `samtools`. Fixtures are
embedded as SAM text and built to BAM with samtools at run time, so the suite is
self-contained — just point it at a samtools binary:

```bash
KIRA_BAM_SAMTOOLS=/path/to/samtools cargo test --release
```

It requires agreement on: `view`/`flagstat`/`idxstats` (byte-identical);
`view -b`, `sort`, and `merge` (identical decoded body). Sort order is a
byte-for-byte contract with samtools, ties included — coordinate sorts on
`(tid, pos, reverse-strand)` with input-stable ties, and `-n` uses natural
(`strnum`) qname order with a READ1/READ2 tie-break. `markdup` must flag the same
duplicate set (via the standard `fixmate -m | sort` prep), `calmd` the same MD/NM
values, and a CRAM round-trip must preserve the mandatory fields. Without
`KIRA_BAM_SAMTOOLS` the golden tests skip, so `cargo test` stays green anywhere.

An optional `KIRA_BAM_FIXTURES=/path/to/bams` folds a directory of real-world
BAMs into the `view`/`flagstat` comparisons.

## License

MIT.
