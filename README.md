# kira-bam

High-performance BAM/SAM/CRAM toolkit written in Rust 2024. Drop-in
samtools-compatible CLI — 31 subcommands, from `view`/`sort`/`index`/`merge`/
`markdup`/`flagstat` through `stats`, `mpileup`, `consensus` and the amplicon
tools ([full list below](#commands)) — plus a library API for embedding directly
into aligners and variant callers, designed to fuse the classic
`bwa | samtools sort | samtools markdup | samtools index` cascade into one
streaming pass.

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

31 subcommands, each mirroring its `samtools` counterpart in name, flag letters
and output format. `kira-bam <command> --help` prints the full flag list.

### Conventions

| Flag | Meaning |
| --- | --- |
| `-@, --threads N` | Worker threads. Default 4 (`sort` defaults to 8). |
| `-T, --reference REF.fa` | Reference FASTA — required to read or write CRAM. The `.fai` must sit next to it. |
| `-o, --output PATH` | Output path. Text-producing commands default to stdout. |
| `-u, --uncompressed` | Write BGZF level 0. |
| `--no-PG` | Do not append the `@PG` provenance record to the output header. |

Input format (SAM / BAM / CRAM) is sniffed from the magic bytes, not from the
file extension — a mismatched extension only prints a warning. Four commands
spell the reference differently, following samtools: `collate` uses `-T` for its
temp directory so its reference is long-only (`--reference`), `mpileup` takes
`-f/--fasta-ref`, `stats` takes `-r/--reference`, and `calmd` takes it as a
positional argument.

### Read, filter, write

```bash
kira-bam view     [-b|-S|-C|-u] [-h|-H] [-f INT] [-F INT] [-q INT] [-c] [-X]
                  [-v REGION]... [-T ref.fa] [-@ N] [-o OUT] [--no-PG] in.bam [region ...]
kira-bam head     [-n INT] [-h INT] [-T ref.fa] [-o OUT] in.bam
kira-bam cat      -o out.bam [-l INT] [-h HEADER] [-T ref.fa] [--no-PG] in1.bam in2.bam ...
kira-bam split    [-f FORMAT] [-d TAG] [-@ N] [-T ref.fa] [--no-PG] in.bam
kira-bam reheader -o out.bam [-T ref.fa] [--no-PG] new_header.sam in.bam
kira-bam addreplacerg -o out.bam -r "ID:foo\tSM:bar" [-m MODE] [-u] [-@ N] [-T ref.fa] [--no-PG] in.bam
kira-bam quickcheck [-v] [-u] in1.bam in2.bam ...
```

- **`view`** — `-b`/`-S`/`-C`/`-u` select the output format (`-C` needs `-T`),
  `-h`/`-H` control the header, `-f`/`-F` take decimal or `0x`-hex flag masks,
  `-q` filters on MAPQ, `-c` counts instead of printing, `-X` drops every aux
  tag, and `-v/--exclude-region` (repeatable) subtracts regions from the output.
  Positional regions need a `.bai`/`.csi` index alongside the input.
- **`head`** — first `-n` records (default 5); `-h` limits the header lines
  printed (0 = all).
- **`cat`** — concatenates BAMs that share a header, without recompressing the
  bodies. `-l -1` keeps the default level, `0` writes uncompressed; `-h` takes
  the header from another file instead of input #0.
- **`split`** — one output file per distinct tag value. `-d` picks the tag
  (default `RG`), `-f` the filename pattern, where `%*` expands to the input
  file stem and `%!` to the tag value (default `%*_%!.bam`).
- **`reheader`** — replaces the header and copies the body blocks untouched.
  `-o` is mandatory: unlike samtools, kira-bam never rewrites in place.
- **`addreplacerg`** — `-r` is the new `@RG` line; `-m` is one of `orphan-only`
  (default), `overwrite_all`, `overwrite_orphans`, `lenient`.
- **`quickcheck`** — exits non-zero if any input is truncated or corrupt. `-v`
  prints one diagnostic per failing file; `-u` treats unrecognized formats as OK
  (samtools fails them).

### Order, merge, deduplicate

```bash
kira-bam sort    [-n] [-u|-l INT] [-m MEM] [--tmpdir DIR] [--reference ref.fa]
                 [-f INT] [-F INT] [--markdup] [-@ N] [-o out.bam] [--no-PG] in.bam
kira-bam collate [-f] [-r BUCKETS] [-T TMPDIR] [--reference ref.fa] [-@ N] [--no-PG] in.bam out_prefix
kira-bam merge   -o out.bam [-n] [-f] [-u] [-T ref.fa] [-@ N] [--no-PG] in1.bam in2.bam ...
kira-bam fixmate [-r] [-m] [-c] [-u] [-@ N] [-T ref.fa] [--no-PG]
                 [--assume-query-grouped] in.namesorted.bam out.bam
kira-bam markdup [-r] [-d INT] [-s STATS] [-S] [-u] [-@ N] [-T ref.fa]
                 [--barcode-tag TAG] [--mode-ancient] [--no-PG] in.sorted.bam out.bam
```

- **`sort`** — coordinate order by default, `-n` for natural qname order. `-m`
  is the per-thread in-memory chunk budget and accepts `768M`, `4G`, a plain
  byte count, or `auto` (3/4 of total RAM, floor 512 MB). `-f`/`-F` drop records
  by flag mask during the pass. `--markdup` folds duplicate marking into the
  sort — same result as piping into `markdup`, one BAM roundtrip cheaper — and
  is tuned by `--markdup-barcode-tag`, `--markdup-mode-ancient` and
  `--markdup-supplementary`; it only applies when the input fits in a single
  in-memory chunk. If the sort would spill, it fails with an error and writes
  no output — raise `-m` (or use `-m auto`), or run `kira-bam markdup`
  separately on the sorted output.
- **`collate`** — groups records of the same qname together without a full
  name-sort. Writes `<out_prefix>.bam`. `-f` collates only primary alignments;
  `-r` sets the hash bucket count (default 1048576 — larger means less spill and
  more RAM); `-T` is the temp directory.
- **`merge`** — merges already-sorted inputs; `-n` for name-sorted inputs, `-f`
  overwrites an existing output.
- **`fixmate`** — repairs mate position, mate ref, MC/MS/MQ tags and the 0x2
  proper-pair bit; input must be queryname-sorted. `-r` drops unmapped and
  secondary/supplementary reads, `-m` adds the MC/MS tags markdup needs, `-c`
  recomputes the CIGAR. `--assume-query-grouped` skips the header sort-order
  check — only correct if mates really are adjacent.
- **`markdup`** — flags FLAG 0x400; `-r` removes the duplicates instead. `-d`
  enables optical-duplicate detection at a pixel distance (cluster coordinates
  are parsed from Illumina read names) and reports the optical subset in the
  `-s` stats file. `--barcode-tag` groups by UMI / molecular index (e.g. `RX`),
  `--mode-ancient` collapses forward and reverse strands together for aDNA, and
  `-S` also flags the supplementary/secondary records of a duplicate template.

### Index

```bash
kira-bam index [-b|-c] [--min-shift N] [--depth N] [-@ N] [-o OUT] in.sorted.bam
kira-bam faidx [--line-width N] [-o OUT] ref.fa [region ...]
kira-bam fqidx [-o OUT] reads.fastq [name ...]
kira-bam dict  [-a ASSEMBLY] [-s SPECIES] [-u URI] [-o OUT] ref.fa
```

- **`index`** — BAI by default, `-c` forces CSI (needed for chromosomes larger
  than 2^29 bp); `--min-shift`/`--depth` tune the CSI binning (defaults 14/5).
- **`faidx`** — builds `ref.fa.fai` if missing, then extracts the requested
  regions (`chr` or `chr:start-end`) as FASTA at `--line-width` columns
  (default 60). With no region it just (re)builds the index.
- **`fqidx`** — the same for FASTQ: indexes the file, and prints each named read
  as a 4-line FASTQ record.
- **`dict`** — writes a `@HD`/`@SQ` sequence dictionary (with M5 checksums) for
  a FASTA, to stdout unless `-o` is given. `-a`/`-s`/`-u` fill in the `AS`, `SP`
  and `UR` fields.

### Statistics and coverage

```bash
kira-bam flagstat      [-O text|json] [-T ref.fa] [-@ N] in.bam
kira-bam idxstats      in.indexed.bam
kira-bam stats         [-r ref.fa] [-t regions.bed] [--sample-name NAME] [-q INT] [-@ N] [-o OUT] in.bam
kira-bam depth         [-a] [-b regions.bed] [-r REGION]... [-q INT] [-Q INT]
                       [-T ref.fa] [-@ N] [-o OUT] in.indexed.bam
kira-bam coverage      [-r REGION]... [-q INT] [-Q INT] [-H] [--columns LIST]
                       [-T ref.fa] [-@ N] [-o OUT] in.bam
kira-bam bedcov        [-Q INT] [-o OUT] regions.bed in1.bam [in2.bam ...]
kira-bam samples       [-X] [-T ref.fa] [-o OUT] in1.bam [in2.bam ...]
kira-bam ampliconstats [-q INT] [-o OUT] amplicons.bed in1.bam [in2.bam ...]
```

- **`flagstat`** — the usual QC-pass/QC-fail counts; `-O json` emits the machine
  form.
- **`idxstats`** — per-reference mapped/unmapped counts read straight from the
  `.bai`/`.csi`, so it needs an indexed BAM and no body scan.
- **`stats`** — the `samtools stats` report. `-r` adds the GC and error-rate
  sections, `-t` restricts to a BED, `-q` sets a MAPQ floor, `--sample-name`
  annotates the output.
- **`depth`** — per-base depth. `-a` also prints zero-coverage positions, `-b`
  restricts to a BED, `-r` (repeatable) to regions, `-q`/`-Q` set MAPQ and base
  quality floors.
- **`coverage`** — one summary row per reference (or per `-r` region). `-H`
  suppresses the header line; `--columns` takes a comma-separated subset of
  `rname,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq`.
- **`bedcov`** — appends one coverage column per input BAM to each BED interval.
- **`samples`** — lists the `SM:` values found in `@RG` headers, one row per
  sample; `-X` prints one row per `@RG` instead.
- **`ampliconstats`** — per-amplicon depth and read counts for the amplicons in
  the BED, across one or more BAMs.

### Pileup and consensus

```bash
kira-bam mpileup   [-f ref.fa] [-r REGION]... [-q INT] [-Q INT] [-d INT] [-a] [-o OUT] in.bam ...
kira-bam consensus [-q INT] [-Q INT] [-d INT] [--call-fraction F] [--line-width N]
                   [-T ref.fa] [-o OUT] in.sorted.bam
kira-bam tview     -p chr:pos [-w INT] [-d INT] [--no-color] in.indexed.bam [ref.fa]
```

- **`mpileup`** — classic pileup text. `-f` supplies the reference column, `-r`
  (repeatable) restricts to regions, `-d` caps depth per position (default
  8000), `-Q` defaults to 13. `-r` bounds every row, covered or zero-filled;
  `-a` emits every position including zero-coverage ones, from the start to the
  end (or the `-r` span) of each reference sequence that has reads — reference
  sequences without reads stay absent even if a region names them. Multiple
  input BAMs each add their own depth/bases/quality triplet, in the order the
  files are given; a position is emitted when any input covers it, and an
  input lacking coverage there contributes a `0 * *` triplet.
- **`consensus`** — majority-call consensus as FASTA. `--call-fraction` is the
  fraction a base must reach to be called (default 0.5, below it becomes `N`),
  `--line-width` the FASTA wrap (default 60). Input must be coordinate-sorted.
- **`tview`** — non-interactive text alignment viewer: prints a fixed-width pile
  around `-p chr:pos` (samtools `tview -d T` style), `-w` bases wide and at most
  `-d` reads deep. `--no-color` suppresses the ANSI escapes.

### Conversion and clipping

```bash
kira-bam fastq        [-1 R1] [-2 R2] [-s SINGLETONS] [-0 ORPHANS] [--fasta] [-N]
                      [--collate] [-T ref.fa] [-@ N] in.bam
kira-bam import       -o out.bam [-s SAMPLE] [-g RG_ID] [-C] [-u] [-@ N] [--no-PG] r1.fq [r2.fq]
kira-bam calmd        [-b] [-u] [-e] [-@ N] [-o OUT] [--no-PG] in.bam ref.fa
kira-bam ampliconclip -b primers.bed -o out.bam [--hard-clip] [-u] [-@ N] [-T ref.fa] [--no-PG] in.bam
```

- **`fastq`** — splits reads into `-1`/`-2`, with `-s` for singletons and `-0`
  for orphans (paired records that are neither READ1 nor READ2). `--fasta` drops
  the quality lines, `-N` appends the CASAVA `/1` `/2` suffixes, `--collate`
  groups mates first when the input is not name-sorted.
- **`import`** — FASTQ to unaligned BAM. `-s` fills `@RG SM`, `-g` the
  read-group ID (default `1`), `-C` parses the CASAVA suffix and pulls the
  barcode into a `BC` tag.
- **`calmd`** — recomputes MD and NM against the reference (a positional
  argument here, indexed with `.fai`). `-b` writes BAM. `-e` is accepted for
  samtools compatibility but is currently a no-op.
- **`ampliconclip`** — soft-clips reads at the primer boundaries in the BED;
  `--hard-clip` hard-clips instead. Input must be coordinate-sorted.

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
