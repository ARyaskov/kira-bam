//! Golden tests: run `kira-bam` and `samtools` on the same inputs and require
//! their outputs to agree.
//!
//! These are gated on a real `samtools` binary, pointed to by the
//! `KIRA_BAM_SAMTOOLS` environment variable. When it is unset every test
//! returns early (prints a skip note), so `cargo test` stays green on machines
//! without samtools. To run them:
//!
//! ```bash
//! KIRA_BAM_SAMTOOLS=/path/to/samtools cargo test --release
//! ```
//!
//! Fixtures are embedded as SAM text and materialised to BAM with samtools at
//! run time, so the suite is fully self-contained. An optional directory of
//! real-world BAMs can additionally be supplied via `KIRA_BAM_FIXTURES`; the
//! `view`/`flagstat` tests fold those in too.
//!
//! ## What "agree" means per command
//! * `view`, `flagstat`, `idxstats` — byte-identical output.
//! * `view -b`, `sort` — identical decoded body (sort also tolerates an
//!   alternative but still-valid coordinate tie-order, verified explicitly).
//! * `markdup` — identical *set* of duplicate-flagged reads (prepared through
//!   the standard `fixmate -m | sort` pipeline both tools require).
//! * `calmd` — identical MD/NM tag values (tag order is allowed to differ).
//! * `mpileup` — CHROM/POS/REF plus every sample's depth/bases/quals triplet
//!   after folding samtools' `^`/`$`/strand boundary markup (kira omits that
//!   markup, and the fixture is deliberately indel-free).
//! * CRAM round-trip — identical mandatory fields (cols 1–11) after decode.
//! * `depth -a -r` — byte-identical output on the shared flag subset (no
//!   `-aa`/`-J`/`-d`/`-s` and no paired/secondary/duplicate reads), plus an
//!   ungated exact-output check on the coordinate-boundary fixture.
//! * `bedcov` — byte-identical echoed BED columns plus one per-input column,
//!   including the missing-`@SQ` row, plus an ungated exact-output check built
//!   entirely with kira (`sort` + `index` feed `bedcov`'s indexed reader).

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::process::Command;

// ── harness plumbing ──────────────────────────────────────────────────────────

fn samtools() -> Option<String> {
    std::env::var("KIRA_BAM_SAMTOOLS")
        .ok()
        .filter(|s| !s.is_empty())
}

fn kira_bam() -> String {
    env!("CARGO_BIN_EXE_kira-bam").to_string()
}

fn fixtures_dir() -> Option<String> {
    std::env::var("KIRA_BAM_FIXTURES")
        .ok()
        .filter(|s| !s.is_empty())
}

fn strip_cr(s: &str) -> String {
    s.replace('\r', "")
}

/// Run a command, asserting success, returning its stdout as a CR-stripped String.
fn run_ok(label: &str, cmd: &mut Command) -> String {
    let out = cmd
        .output()
        .unwrap_or_else(|e| panic!("{label}: spawn failed: {e}"));
    assert!(
        out.status.success(),
        "{label}: exit {:?}\nstderr: {}",
        out.status.code(),
        String::from_utf8_lossy(&out.stderr)
    );
    strip_cr(&String::from_utf8_lossy(&out.stdout))
}

/// Materialise embedded SAM text into a BAM via samtools (the neutral builder).
fn build_bam(samtools: &str, dir: &Path, name: &str, sam: &str) -> PathBuf {
    let sam_path = dir.join(format!("{name}.sam"));
    std::fs::write(&sam_path, sam).unwrap();
    let bam_path = dir.join(format!("{name}.bam"));
    run_ok(
        "build_bam",
        Command::new(samtools)
            .args(["view", "-b", "-o"])
            .arg(&bam_path)
            .arg(&sam_path),
    );
    bam_path
}

/// Decode a BAM/CRAM body (no header) through samtools, so the *decoder* is held
/// constant and only the encoded content is compared.
fn decode(samtools: &str, bam: &Path, extra: &[&str]) -> String {
    run_ok(
        "decode",
        Command::new(samtools).arg("view").args(extra).arg(bam),
    )
}

fn tmpdir() -> tempfile::TempDir {
    tempfile::TempDir::new().unwrap()
}

// ── embedded fixtures ─────────────────────────────────────────────────────────

/// Variety pack: two references, proper pairs, soft/hard clips, an indel, a
/// deletion, secondary + supplementary records, and an unmapped pair.
const BASIC_SAM: &str = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:2000
@SQ\tSN:chr2\tLN:2000
rA\t99\tchr1\t100\t60\t10M\t=\t300\t210\tACGTACGTAC\tIIIIIIIIII
rA\t147\tchr1\t300\t60\t8M2S\t=\t100\t-210\tACGTACGTAC\tFFFFFFFFFF
rB\t99\tchr1\t120\t30\t2S8M\t=\t350\t240\tTTTTTTTTTT\tIIIIIIIIII
rB\t147\tchr1\t350\t30\t10M\t=\t120\t-240\tTTTTTTTTTT\tIIIIIIIIII
rC\t0\tchr1\t150\t60\t5M2D5M\t*\t0\t0\tAACCGGTTAA\tIIIIIIIIII
rE\t256\tchr1\t500\t0\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
rF\t2048\tchr1\t600\t60\t5H10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
rH\t0\tchr2\t200\t60\t4M1I5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
rI\t77\t*\t0\t0\t*\t*\t0\t0\tNNNNNNNNNN\tIIIIIIIIII
rJ\t141\t*\t0\t0\t*\t*\t0\t0\tNNNNNNNNNN\tIIIIIIIIII
";

/// Name-collated read pairs where f2 is a PCR duplicate of f1 (identical
/// unclipped coordinates) and f3 is distinct. Soft clips ensure the unclipped-5′
/// logic is exercised.
const PAIRS_SAM: &str = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:2000
f1\t99\tchr1\t100\t60\t2S8M\t=\t300\t210\tACGTACGTAC\tIIIIIIIIII
f1\t147\tchr1\t300\t60\t8M2S\t=\t100\t-210\tACGTACGTAC\tIIIIIIIIII
f2\t99\tchr1\t100\t60\t2S8M\t=\t300\t210\tACGTACGTAC\tFFFFFFFFFF
f2\t147\tchr1\t300\t60\t8M2S\t=\t100\t-210\tACGTACGTAC\tFFFFFFFFFF
f3\t99\tchr1\t500\t60\t10M\t=\t700\t210\tACGTACGTAC\tIIIIIIIIII
f3\t147\tchr1\t700\t60\t10M\t=\t500\t-210\tACGTACGTAC\tIIIIIIIIII
";

/// Single-end reads aligned to `reference_fasta()` with a mismatch, a deletion,
/// and an insertion — exercises MD/NM and CRAM reference encoding.
const ALIGNED_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:200
m1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGGACGTAC\tIIIIIIIIII
m2\t0\tchr1\t21\t60\t5M2D5M\t*\t0\t0\tACGTAACGTA\tIIIIIIIIII
m3\t0\tchr1\t41\t60\t4M1I5M\t*\t0\t0\tACGTTCGTAC\tIIIIIIIIII
";

/// Many records sharing (tid, pos) in scrambled order, with both strands — pins
/// the coordinate tie-break (forward before reverse, then input order).
const TIES_SAM: &str = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:2000
zebra\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
alpha\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
mike\t16\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
alpha\t16\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
bravo\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
delta\t0\tchr1\t200\t60\t5M\t*\t0\t0\tACGTA\tIIIII
";

/// Numeric-tail qnames out of order, plus one name shared by read1/read2/supp —
/// pins natural (`strnum`) qname order and the READ1/READ2 tie-break.
const NAME_TIES_SAM: &str = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:2000
r10\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
r2\t147\tchr1\t500\t60\t5M\t=\t100\t-405\tACGTA\tIIIII
r2\t99\tchr1\t100\t60\t5M\t=\t500\t405\tACGTA\tIIIII
r2\t2048\tchr1\t300\t60\t5M\t*\t0\t0\tACGTA\tIIIII
r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
";

/// Two coordinate-sorted inputs for `merge`, with records at identical positions
/// across sources (cross-source ties) and both strands.
const MERGE_A_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
a1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
a2\t0\tchr1\t200\t60\t5M\t*\t0\t0\tACGTA\tIIIII
";

const MERGE_B_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
b1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
b2\t16\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
b3\t0\tchr1\t200\t60\t5M\t*\t0\t0\tACGTA\tIIIII
";

/// 60 bp `ACGT`-repeat FASTA matching the mpileup fixtures (one contig, so a
/// single `faidx` write enables both tools' reference styling).
fn mpileup_reference_fasta() -> String {
    let seq: String = (0..60).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect();
    format!(">chr1\n{seq}\n")
}

/// mpileup golden, input A: one forward 6M read at pos 1, all matching the
/// `ACGT` repeat, quals `IJKLMN` (Phred 40–45, above the shared `-Q 13` floor).
const MPILEUP_A_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:60
a1\t0\tchr1\t1\t60\t6M\t*\t0\t0\tACGTAC\tIJKLMN
";

/// mpileup golden, input B: b1 covers pos 2–7 (`CGTACG`, all matching, quals
/// `ABCDEF`); b2 covers pos 6–11 (`CGTTCG`, one mismatch at pos 9: reads `T`
/// against reference `A`, quals `OPQRST`). Forward-only, indel-free, so the
/// only expected samtools markup is read-start `^`+mapQ and read-end `$`.
const MPILEUP_B_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:60
b1\t0\tchr1\t2\t60\t6M\t*\t0\t0\tCGTACG\tABCDEF
b2\t0\tchr1\t6\t60\t6M\t*\t0\t0\tCGTTCG\tOPQRST
";

/// 200 bp `ACGT`-repeat reference matching `ALIGNED_SAM`.
fn reference_fasta() -> String {
    let seq: String = (0..200).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect();
    format!(">chr1\n{seq}\n")
}

/// All BAMs to feed the format/aggregate tests: embedded BASIC plus any
/// real-world fixtures supplied via KIRA_BAM_FIXTURES.
fn corpus(samtools: &str, dir: &Path) -> Vec<PathBuf> {
    let mut bams = vec![build_bam(samtools, dir, "basic", BASIC_SAM)];
    if let Some(fix) = fixtures_dir() {
        bams.extend(list_bams(Path::new(&fix)));
    }
    bams
}

fn list_bams(dir: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    if let Ok(entries) = std::fs::read_dir(dir) {
        for e in entries.flatten() {
            let p = e.path();
            if p.is_dir() {
                out.extend(list_bams(&p));
            } else if p.extension().and_then(|s| s.to_str()) == Some("bam") {
                out.push(p);
            }
        }
    }
    out
}

macro_rules! require_samtools {
    () => {
        match samtools() {
            Some(s) => s,
            None => {
                eprintln!("skip: set KIRA_BAM_SAMTOOLS to run golden tests");
                return;
            }
        }
    };
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[test]
fn flagstat_byte_identical() {
    let st = require_samtools!();
    let dir = tmpdir();
    for bam in corpus(&st, dir.path()) {
        let s = run_ok(
            "samtools flagstat",
            Command::new(&st).arg("flagstat").arg(&bam),
        );
        let k = run_ok(
            "kira flagstat",
            Command::new(kira_bam()).arg("flagstat").arg(&bam),
        );
        assert_eq!(s, k, "flagstat mismatch for {}", bam.display());
    }
}

#[test]
fn view_body_byte_identical() {
    let st = require_samtools!();
    let dir = tmpdir();
    for bam in corpus(&st, dir.path()) {
        let s = run_ok("samtools view", Command::new(&st).arg("view").arg(&bam));
        let k = run_ok("kira view", Command::new(kira_bam()).arg("view").arg(&bam));
        assert_eq!(s, k, "view body mismatch for {}", bam.display());
    }
}

#[test]
fn view_bam_roundtrip_identical() {
    let st = require_samtools!();
    let dir = tmpdir();
    let bam = build_bam(&st, dir.path(), "basic", BASIC_SAM);
    // kira encodes BAM; decode with samtools and compare to the original body.
    let out = dir.path().join("kira.bam");
    run_ok(
        "kira view -b",
        Command::new(kira_bam())
            .args(["view", "-b", "-o"])
            .arg(&out)
            .arg(&bam),
    );
    assert_eq!(
        decode(&st, &bam, &[]),
        decode(&st, &out, &[]),
        "kira BAM round-trip changed the body"
    );
}

#[test]
fn sort_matches_samtools() {
    let st = require_samtools!();
    let dir = tmpdir();
    let bam = build_bam(&st, dir.path(), "basic", BASIC_SAM);

    let ks = dir.path().join("k.sorted.bam");
    run_ok(
        "kira sort",
        Command::new(kira_bam())
            .args(["sort", "-o"])
            .arg(&ks)
            .arg(&bam),
    );
    let ss = dir.path().join("s.sorted.bam");
    run_ok(
        "samtools sort",
        Command::new(&st).args(["sort", "-o"]).arg(&ss).arg(&bam),
    );

    let kbody = decode(&st, &ks, &[]);
    let sbody = decode(&st, &ss, &[]);
    if kbody == sbody {
        return; // exact tie-order match — ideal
    }
    // Otherwise the records must still be the same multiset AND kira's output a
    // valid coordinate sort (only the order of equal-coordinate records differs).
    assert_eq!(
        multiset(&kbody),
        multiset(&sbody),
        "sort changed the set of records"
    );
    assert_coordinate_sorted(&st, &ks);
}

#[test]
fn sort_coordinate_ties_byte_identical() {
    let st = require_samtools!();
    let dir = tmpdir();
    let bam = build_bam(&st, dir.path(), "ties", TIES_SAM);
    let ks = dir.path().join("k.bam");
    run_ok(
        "kira sort",
        Command::new(kira_bam())
            .args(["sort", "-o"])
            .arg(&ks)
            .arg(&bam),
    );
    let ss = dir.path().join("s.bam");
    run_ok(
        "samtools sort",
        Command::new(&st).args(["sort", "-o"]).arg(&ss).arg(&bam),
    );
    // Strict: with the htslib-matching (tid, pos, reverse) key + stable order, the
    // coordinate tie-break must be byte-identical, not merely the same multiset.
    assert_eq!(
        decode(&st, &ss, &[]),
        decode(&st, &ks, &[]),
        "coordinate tie-order differs from samtools"
    );
}

#[test]
fn sort_name_order_byte_identical() {
    let st = require_samtools!();
    let dir = tmpdir();
    let bam = build_bam(&st, dir.path(), "nties", NAME_TIES_SAM);
    let ks = dir.path().join("k.bam");
    run_ok(
        "kira sort -n",
        Command::new(kira_bam())
            .args(["sort", "-n", "-o"])
            .arg(&ks)
            .arg(&bam),
    );
    let ss = dir.path().join("s.bam");
    run_ok(
        "samtools sort -n",
        Command::new(&st)
            .args(["sort", "-n", "-o"])
            .arg(&ss)
            .arg(&bam),
    );
    assert_eq!(
        decode(&st, &ss, &[]),
        decode(&st, &ks, &[]),
        "name sort order differs from samtools (strnum / READ1-READ2 tie-break)"
    );
}

#[test]
fn merge_byte_identical() {
    let st = require_samtools!();
    let dir = tmpdir();
    let a = build_bam(&st, dir.path(), "ma", MERGE_A_SAM);
    let b = build_bam(&st, dir.path(), "mb", MERGE_B_SAM);
    // Inputs are already coordinate-sorted; merge with both and compare.
    let km = dir.path().join("k.merged.bam");
    run_ok(
        "kira merge",
        Command::new(kira_bam())
            .args(["merge", "-f", "-o"])
            .arg(&km)
            .arg(&a)
            .arg(&b),
    );
    let sm = dir.path().join("s.merged.bam");
    run_ok(
        "samtools merge",
        Command::new(&st)
            .arg("merge")
            .arg("-f")
            .arg(&sm)
            .arg(&a)
            .arg(&b),
    );
    assert_eq!(
        decode(&st, &sm, &[]),
        decode(&st, &km, &[]),
        "merge produced a different order than samtools"
    );
}

#[test]
fn index_and_idxstats_match() {
    let st = require_samtools!();
    let dir = tmpdir();
    let bam = build_bam(&st, dir.path(), "basic", BASIC_SAM);
    let sorted = dir.path().join("sorted.bam");
    run_ok(
        "samtools sort",
        Command::new(&st)
            .args(["sort", "-o"])
            .arg(&sorted)
            .arg(&bam),
    );

    // Build the index with kira, then make samtools read it: tests both the
    // kira-written .bai and kira idxstats against samtools.
    run_ok(
        "kira index",
        Command::new(kira_bam()).arg("index").arg(&sorted),
    );
    let s = run_ok(
        "samtools idxstats",
        Command::new(&st).arg("idxstats").arg(&sorted),
    );
    let k = run_ok(
        "kira idxstats",
        Command::new(kira_bam()).arg("idxstats").arg(&sorted),
    );
    assert_eq!(s, k, "idxstats mismatch on kira-built index");
}

#[test]
fn markdup_duplicate_set_matches() {
    let st = require_samtools!();
    let dir = tmpdir();
    let raw = build_bam(&st, dir.path(), "pairs", PAIRS_SAM);

    // Both tools require the standard prep: fixmate -m (adds MC/ms) then sort.
    // Prepare once with samtools so the markdup *input* is identical for both.
    let fm = dir.path().join("fm.bam");
    run_ok(
        "samtools fixmate",
        Command::new(&st).args(["fixmate", "-m"]).arg(&raw).arg(&fm),
    );
    let prepared = dir.path().join("prepared.bam");
    run_ok(
        "samtools sort",
        Command::new(&st)
            .args(["sort", "-o"])
            .arg(&prepared)
            .arg(&fm),
    );

    let s_md = dir.path().join("s.md.bam");
    run_ok(
        "samtools markdup",
        Command::new(&st).arg("markdup").arg(&prepared).arg(&s_md),
    );
    let k_md = dir.path().join("k.md.bam");
    run_ok(
        "kira markdup",
        Command::new(kira_bam())
            .arg("markdup")
            .arg(&prepared)
            .arg(&k_md),
    );

    assert_eq!(
        dup_flagged_qnames(&st, &s_md),
        dup_flagged_qnames(&st, &k_md),
        "markdup flagged a different set of duplicates"
    );
}

#[test]
fn calmd_md_nm_values_match() {
    let st = require_samtools!();
    let dir = tmpdir();
    let refp = dir.path().join("ref.fa");
    std::fs::write(&refp, reference_fasta()).unwrap();
    run_ok("samtools faidx", Command::new(&st).arg("faidx").arg(&refp));
    let bam = build_bam(&st, dir.path(), "aln", ALIGNED_SAM);

    let s = run_ok(
        "samtools calmd",
        Command::new(&st).arg("calmd").arg(&bam).arg(&refp),
    );
    let k = run_ok(
        "kira calmd",
        Command::new(kira_bam()).arg("calmd").arg(&bam).arg(&refp),
    );
    // Compare MD/NM by read, order-independent (samtools emits NM,MD; kira MD,NM).
    assert_eq!(
        md_nm_by_read(&s),
        md_nm_by_read(&k),
        "calmd MD/NM values differ"
    );
}

#[test]
fn cram_roundtrip_preserves_alignment() {
    let st = require_samtools!();
    let dir = tmpdir();
    let refp = dir.path().join("ref.fa");
    std::fs::write(&refp, reference_fasta()).unwrap();
    run_ok("samtools faidx", Command::new(&st).arg("faidx").arg(&refp));
    let bam = build_bam(&st, dir.path(), "aln", ALIGNED_SAM);

    let cram = dir.path().join("k.cram");
    run_ok(
        "kira view -C",
        Command::new(kira_bam())
            .args(["view", "-C", "-T"])
            .arg(&refp)
            .args(["-o"])
            .arg(&cram)
            .arg(&bam),
    );
    // CRAM is reference-encoded; samtools regenerates MD/NM on decode, so compare
    // only the mandatory fields (cols 1–11), which must round-trip losslessly.
    let original = mandatory_fields(&decode(&st, &bam, &[]));
    let restored = mandatory_fields(&decode(&st, &cram, &["-T", refp.to_str().unwrap()]));
    assert_eq!(
        original, restored,
        "CRAM round-trip changed an alignment field"
    );
}

/// Depth edge fixture: one contig (`chr1` LN:8) with single-base reads at the
/// first and last 1-based coordinate positions, CIGAR 1M. Single-end, flags 0,
/// MAPQ 60, base quality `I` — above all defaults (excl 0x704, min MAPQ 0,
/// min BQ 0) and free of D/N and overlapping pairs, so every observation
/// counts exactly once.
const DEPTH_EDGE_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:8
d1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI
d8\t0\tchr1\t8\t60\t1M\t*\t0\t0\tA\tI
";

/// The exact `depth -a -r chr1:1-8` output both tools must produce for
/// [`DEPTH_EDGE_SAM`]: every contig position 1–8 (1-based inclusive), depth 1
/// at the read positions and 0 elsewhere.
fn depth_edge_expected() -> String {
    (1..=8)
        .map(|p| format!("chr1\t{p}\t{}", u32::from(p == 1 || p == 8)))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

#[test]
fn depth_region_all_positions_exact() {
    // Ungated: kira reads plain SAM text directly, so the first/last-coordinate
    // boundary behavior is checked even without samtools.
    let dir = tmpdir();
    let sam = dir.path().join("depth_edge.sam");
    std::fs::write(&sam, DEPTH_EDGE_SAM).unwrap();
    let out = run_ok(
        "kira depth -a -r",
        Command::new(kira_bam())
            .args(["depth", "-a", "-r", "chr1:1-8"])
            .arg(&sam),
    );
    assert_eq!(out, depth_edge_expected(), "kira depth edge mismatch");
}

#[test]
fn depth_region_all_positions_matches_samtools() {
    let st = require_samtools!();
    let dir = tmpdir();
    let bam = build_bam(&st, dir.path(), "depth_edge", DEPTH_EDGE_SAM);
    // `-r` region access goes through the index.
    run_ok("samtools index", Command::new(&st).arg("index").arg(&bam));
    let s = run_ok(
        "samtools depth -a -r",
        Command::new(&st)
            .args(["depth", "-a", "-r", "chr1:1-8"])
            .arg(&bam),
    );
    let k = run_ok(
        "kira depth -a -r",
        Command::new(kira_bam())
            .args(["depth", "-a", "-r", "chr1:1-8"])
            .arg(&bam),
    );
    assert_eq!(s, k, "depth -a -r differs from samtools");
    assert_eq!(k, depth_edge_expected());
}

// ── bedcov ────────────────────────────────────────────────────────────────────

/// `bedcov` fixture: one 20 bp contig (`chr1`) and three single-end all-`M`
/// reads in coordinate order, flag 0, MAPQ 60, qual `I` — above every default
/// filter, with no D/N, MAPQ 255, pairs, secondary or duplicate records, so
/// each read contributes exactly its covered reference span. The recordings
/// deliberately leave `chr2` out of the header so the missing-`@SQ` path (0 +
/// a stderr warning) is exercised.
///
/// Covered spans (BED 0-based half-open): b1 [2,6), b2 [4,8), b3 [14,20).
const BEDCOV_SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:20
b1\t0\tchr1\t3\t60\t4M\t*\t0\t0\tACGT\tIIII
b2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII
b3\t0\tchr1\t15\t60\t6M\t*\t0\t0\tACGTAC\tIIIIII
";

/// BED rows compared against samtools: the four `chr1` intervals only. The
/// missing-`@SQ` `chr2` row is covered only by the ungated Kira expected-
/// output test `bedcov_rows_exact_ungated`, because cross-tool behavior for
/// a BED name absent from `@SQ` isn't established (samtools not pinned);
/// `bedcov_rows_match_samtools` therefore never sees that row.
const BEDCOV_BED_SHARED: &str = "chr1\t0\t10\n\
                                 chr1\t4\t5\n\
                                 chr1\t6\t10\n\
                                 chr1\t14\t20\n";

/// Tab-separated BED used with [`BEDCOV_SAM`] in the ungated Kira comparison.
/// Coordinates are 0-based half-open; the trailing Kira-only row
/// `chr2 0 10` names a contig no input has in its header.
/// Built at runtime: `concat!` accepts literals, not const identifiers.
fn bedcov_bed() -> String {
    BEDCOV_BED_SHARED.to_owned() + "chr2\t0\t10\n"
}

/// Expected `bedcov` output for [`BEDCOV_BED_SHARED`] + [`BEDCOV_SAM`]: the BED
/// columns echoed verbatim plus one sum-of-covered-base column per input.
///
/// `chr1 0 10` (1-based 1–10): b1 contributes 4 bases (3–6) and b2 contributes
/// 4 (5–8); their overlap 5–6 counts twice → 8. `chr1 4 5` is the single base 5 under both b1 and b2
/// → 2. `chr1 6 10` (1-based 7–10) reaches only b2's tail (7–8) — b1's last
/// base 6 is excluded by the half-open start → 2. `chr1 14 20` is 1-based
/// 15–20, entirely under b3 → 6.
const BEDCOV_EXPECTED_SHARED: &str = "chr1\t0\t10\t8\n\
                                      chr1\t4\t5\t2\n\
                                      chr1\t6\t10\t2\n\
                                      chr1\t14\t20\t6\n";

/// `bedcov` output for the ungated Kira comparison: [`BEDCOV_BED_SHARED`] +
/// [`BEDCOV_SAM`] must produce exactly
/// [`BEDCOV_EXPECTED_SHARED`], plus the Kira-only `chr2` row with no
/// `@SQ` anywhere → 0.
fn bedcov_expected() -> String {
    BEDCOV_EXPECTED_SHARED.to_owned() + "chr2\t0\t10\t0\n"
}

#[test]
fn bedcov_rows_exact_ungated() {
    // Ungated: the fixture is plain SAM text, and the whole chain (sort →
    // index → bedcov) runs on kira itself — no samtools needed.
    let dir = tmpdir();
    let sam = dir.path().join("bedcov.sam");
    std::fs::write(&sam, BEDCOV_SAM).unwrap();
    let bed = dir.path().join("bedcov.bed");
    std::fs::write(&bed, bedcov_bed()).unwrap();

    let sorted = dir.path().join("bedcov.sorted.bam");
    run_ok(
        "kira sort",
        Command::new(kira_bam())
            .args(["sort", "-o"])
            .arg(&sorted)
            .arg(&sam),
    );
    run_ok(
        "kira index",
        Command::new(kira_bam()).arg("index").arg(&sorted),
    );
    let out = run_ok(
        "kira bedcov",
        Command::new(kira_bam())
            .arg("bedcov")
            .arg(&bed)
            .arg(&sorted),
    );
    assert_eq!(out, bedcov_expected(), "kira bedcov rows differ");
}

#[test]
fn bedcov_rows_match_samtools() {
    let st = require_samtools!();
    let dir = tmpdir();
    // Shared rows only; the missing-`@SQ` row stays in the Kira-only
    // test above (cross-tool behavior not established).
    let bed = dir.path().join("bedcov.bed");
    std::fs::write(&bed, BEDCOV_BED_SHARED).unwrap();
    let bam = build_bam(&st, dir.path(), "bedcov", BEDCOV_SAM);
    // Both tools take the indexed-BAM route.
    run_ok("samtools index", Command::new(&st).arg("index").arg(&bam));
    let s = run_ok(
        "samtools bedcov",
        Command::new(&st).arg("bedcov").arg(&bed).arg(&bam),
    );
    let k = run_ok(
        "kira bedcov",
        Command::new(kira_bam()).arg("bedcov").arg(&bed).arg(&bam),
    );
    assert_eq!(s, k, "bedcov differs from samtools");
    assert_eq!(
        k, BEDCOV_EXPECTED_SHARED,
        "bedcov rows differ from expected"
    );
}

#[test]
fn mpileup_multi_input_columns_match_samtools() {
    let st = require_samtools!();
    let dir = tmpdir();
    let refp = dir.path().join("ref.fa");
    std::fs::write(&refp, mpileup_reference_fasta()).unwrap();
    run_ok("samtools faidx", Command::new(&st).arg("faidx").arg(&refp));
    let a = build_bam(&st, dir.path(), "mpa", MPILEUP_A_SAM);
    let b = build_bam(&st, dir.path(), "mpb", MPILEUP_B_SAM);

    let k_out = dir.path().join("k.mpileup.tsv");
    run_ok(
        "kira mpileup",
        Command::new(kira_bam())
            .args(["mpileup", "-f"])
            .arg(&refp)
            .args(["-o"])
            .arg(&k_out)
            .args(["-d", "8000", "-q", "0", "-Q", "13"])
            .arg(&a)
            .arg(&b),
    );
    let kira = mpileup_column_projection(&strip_cr(&std::fs::read_to_string(&k_out).unwrap()))
        .unwrap_or_else(|e| panic!("kira mpileup output rejected by the normalizer: {e}"));

    // -B drops BAQ (kira has none) and the filters are matched, so the only
    // remaining representational difference is samtools' boundary/strand
    // markup, which the projection folds away per sample.
    let s_out = run_ok(
        "samtools mpileup",
        Command::new(&st)
            .args(["mpileup", "-B", "-f"])
            .arg(&refp)
            .args(["-d", "8000", "-q", "0", "-Q", "13"])
            .arg(&a)
            .arg(&b),
    );
    let sam = mpileup_column_projection(&s_out)
        .unwrap_or_else(|e| panic!("samtools mpileup output rejected by the normalizer: {e}"));

    // Full per-sample comparison — CHROM/POS/REF plus each input's depth,
    // normalized base calls, and quality bytes in CLI sample order — not just
    // the depth projection.
    assert_eq!(
        sam, kira,
        "mpileup multi-input columns differ from samtools"
    );
}

// ── comparison helpers ────────────────────────────────────────────────────────

/// Parse one sample's `depth bases quals` triplet into normalized
/// `depth bases quals`, or a description of what could not be represented.
///
/// Samtools can emit markup kira omits: the read-start `^` followed by exactly
/// one mapping-quality byte, the read-end `$`, and strand notation (`,`, and
/// lowercase letters, for reverse-strand observations). Those marker bytes take
/// no base or quality slot, so they are folded away while the retained bases
/// keep their pairing order with the quality bytes. Zero coverage is the exact
/// `0 * *` sentinel. Indel (`+`/`-`) and deletion/ref-skip (`*#<>`) tokens are
/// rejected loudly: the golden fixture is deliberately indel-free, so a hit
/// means counting characters would be wrong (one observation must pair with
/// exactly one quality byte), not silently stripped.
fn mpileup_sample_projection(
    depth: &str,
    bases: &str,
    quals: &str,
) -> Result<(String, String, String), String> {
    let depth: u32 = depth
        .parse()
        .map_err(|_| format!("non-numeric depth field: {depth:?}"))?;
    if depth == 0 {
        if bases != "*" || quals != "*" {
            return Err(format!(
                "zero-depth sample is not the `0 * *` sentinel: {bases:?} {quals:?}"
            ));
        }
        return Ok((depth.to_string(), bases.to_string(), quals.to_string()));
    }

    let b = bases.as_bytes();
    let mut normalized = String::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'^' => {
                // Read-start marker: `^` plus exactly one mapQ byte, neither
                // of which is a base observation.
                if i + 1 >= b.len() {
                    return Err("dangling read-start `^` marker".to_string());
                }
                i += 2;
            }
            b'$' => i += 1, // read-end marker, not a base observation
            b'.' | b',' => {
                normalized.push('.');
                i += 1;
            }
            c if c.is_ascii_alphabetic() => {
                normalized.push(c.to_ascii_uppercase() as char);
                i += 1;
            }
            b'+' | b'-' => {
                return Err(format!(
                    "indel markup `{}` in the bases column: the mpileup golden \
                     fixture is indel-free, so this row cannot be paired with \
                     its qualities",
                    b[i] as char
                ));
            }
            b'*' | b'#' | b'<' | b'>' => {
                return Err(format!(
                    "deletion/reference-skip markup `{}` in the bases column: \
                     the mpileup golden fixture has no skipped observations",
                    b[i] as char
                ));
            }
            c => {
                return Err(format!("unsupported pileup base token {:#04x}", c));
            }
        }
    }

    // One base observation must pair with exactly one quality byte, and the
    // depth column must agree with both.
    let n_quals = quals.len();
    if normalized.len() != n_quals || depth as usize != n_quals {
        return Err(format!(
            "bases, qualities and depth are not one-to-one: {} normalized \
             bases, {} quality bytes, depth {depth}",
            normalized.len(),
            n_quals
        ));
    }
    Ok((depth.to_string(), normalized, quals.to_string()))
}

/// CHROM/POS/REF plus every sample's triplet normalized through
/// [`mpileup_sample_projection`], in CLI sample order — the multi-input
/// invariant both tools must keep (depth, base calls and qualities per input).
fn mpileup_column_projection(out: &str) -> Result<Vec<String>, String> {
    out.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            // CHROM/POS/REF plus one depth/bases/quals triplet per input:
            // exactly 3 + 3N fields.
            if f.len() < 6 || !f.len().is_multiple_of(3) {
                return Err(format!("malformed mpileup row: {l:?}"));
            }
            let mut parts: Vec<String> = f[..3].iter().map(|s| (*s).to_owned()).collect();
            for s in f[3..].as_chunks::<3>().0 {
                let (d, b, q) = mpileup_sample_projection(s[0], s[1], s[2])?;
                parts.push(d);
                parts.push(b);
                parts.push(q);
            }
            Ok(parts.join("\t"))
        })
        .collect()
}

#[test]
fn mpileup_normalizer_strips_read_start_and_end_markers() {
    // `^` consumes exactly the next byte (mapQ 60 → `]`); `$` trails the read's
    // final base and takes no base or quality slot.
    let (d, b, q) = mpileup_sample_projection("1", "^].", "I").unwrap();
    assert_eq!((d.as_str(), b.as_str(), q.as_str()), ("1", ".", "I"));
    // Two observations in one row where the first read ends (`.$. ` shape:
    // b1's base, its `$`, then the next read's base).
    let (d, b, q) = mpileup_sample_projection("2", "^].$.", "IJ").unwrap();
    assert_eq!((d.as_str(), b.as_str(), q.as_str()), ("2", "..", "IJ"));
}

#[test]
fn mpileup_normalizer_normalizes_strand_notation() {
    // samtools marks reverse-strand matches with `,` and mismatches lowercase;
    // kira only ever emits `.` and uppercase letters.
    let (d, b, q) = mpileup_sample_projection("5", ",ACac", "!!!!!").unwrap();
    assert_eq!(
        (d.as_str(), b.as_str(), q.as_str()),
        ("5", ".ACAC", "!!!!!")
    );
}

#[test]
fn mpileup_normalizer_requires_the_zero_sentinel() {
    let (d, b, q) = mpileup_sample_projection("0", "*", "*").unwrap();
    assert_eq!((d.as_str(), b.as_str(), q.as_str()), ("0", "*", "*"));
    assert!(mpileup_sample_projection("0", ".", "I").is_err());
}

#[test]
fn mpileup_normalizer_rejects_indel_and_refskip_markup_and_bad_counts() {
    // Indel and deletion/ref-skip tokens must fail clearly, not be stripped.
    for bad in ['+', '-', '*', '#', '<', '>'] {
        let err = mpileup_sample_projection("1", &bad.to_string(), "I").unwrap_err();
        assert!(!err.is_empty(), "no rejection for {bad:?}");
    }
    // A `^` without its mapQ byte is malformed.
    assert!(mpileup_sample_projection("1", "^", "I").is_err());
    // Bases, qualities and depth must be one-to-one.
    assert!(mpileup_sample_projection("2", "..", "I").is_err());
    assert!(mpileup_sample_projection("1", "..", "II").is_err());
    assert!(mpileup_sample_projection("3", "..", "II").is_err());
}

#[test]
fn mpileup_projection_normalizes_samtools_markup_to_kira_shape() {
    // Row for pos 7 of the mpileup golden fixture: kira writes `..` where
    // samtools writes `.$. ` (b's first read ends mid-column), and input A
    // contributes the zero sentinel. Both must project identically.
    let kira = "chr1\t7\tG\t2\t..\tFP\t0\t*\t*";
    let samtools = "chr1\t7\tG\t2\t.$.\tFP\t0\t*\t*";
    let k = mpileup_column_projection(kira).unwrap_or_else(|e| panic!("{e}"));
    let s = mpileup_column_projection(samtools).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(s, k);
    // Exactly 3 + 3N fields: a 1-sample row is accepted, anything else malformed.
    assert!(mpileup_column_projection("chr1\t7\tG\t1\t.\tI").is_ok());
    assert!(mpileup_column_projection("chr1\t7\tG\t2\t..").is_err());
    assert!(mpileup_column_projection("chr1\t7\tG\t2\t..\tFP\t0").is_err());
}

/// Sorted multiset of non-empty lines (order-insensitive record-set equality).
fn multiset(body: &str) -> Vec<String> {
    let mut v: Vec<String> = body
        .lines()
        .filter(|l| !l.is_empty())
        .map(str::to_string)
        .collect();
    v.sort();
    v
}

/// First 11 (mandatory SAM) fields of each record line, dropping optional tags.
fn mandatory_fields(body: &str) -> Vec<String> {
    body.lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').take(11).collect::<Vec<_>>().join("\t"))
        .collect()
}

/// QNAMEs of records carrying the duplicate flag (0x400), as a set.
fn dup_flagged_qnames(samtools: &str, bam: &Path) -> BTreeSet<String> {
    run_ok(
        "samtools view -f1024",
        Command::new(samtools).args(["view", "-f", "1024"]).arg(bam),
    )
    .lines()
    .filter(|l| !l.is_empty())
    .map(|l| l.split('\t').next().unwrap_or("").to_string())
    .collect()
}

/// Map QNAME → (MD, NM) by scanning optional fields, ignoring tag order.
fn md_nm_by_read(body: &str) -> BTreeMap<String, (String, String)> {
    let mut m = BTreeMap::new();
    for line in body.lines() {
        if line.is_empty() || line.starts_with('@') {
            continue;
        }
        let qname = line.split('\t').next().unwrap_or("").to_string();
        let mut md = String::new();
        let mut nm = String::new();
        for f in line.split('\t') {
            if let Some(v) = f.strip_prefix("MD:Z:") {
                md = v.to_string();
            } else if let Some(v) = f.strip_prefix("NM:i:") {
                nm = v.to_string();
            }
        }
        m.insert(qname, (md, nm));
    }
    m
}

/// Assert a BAM is in coordinate order: (reference index, position) is
/// non-decreasing, with unmapped (`*`) records placed last.
fn assert_coordinate_sorted(samtools: &str, bam: &Path) {
    let header = run_ok(
        "samtools view -H",
        Command::new(samtools).args(["view", "-H"]).arg(bam),
    );
    let mut order: HashMap<String, usize> = HashMap::new();
    for line in header.lines() {
        if let Some(rest) = line.strip_prefix("@SQ\t")
            && let Some(sn) = rest.split('\t').find_map(|f| f.strip_prefix("SN:"))
        {
            let idx = order.len();
            order.insert(sn.to_string(), idx);
        }
    }
    let body = decode(samtools, bam, &[]);
    let mut last = (0usize, 0i64);
    for line in body.lines().filter(|l| !l.is_empty()) {
        let mut f = line.split('\t');
        let rname = f.nth(2).unwrap_or("*");
        let pos: i64 = f.next().unwrap_or("0").parse().unwrap_or(0);
        let tid = if rname == "*" {
            usize::MAX
        } else {
            *order.get(rname).unwrap_or(&usize::MAX)
        };
        let cur = (tid, pos);
        assert!(cur >= last, "not coordinate-sorted at {line}");
        last = cur;
    }
}
