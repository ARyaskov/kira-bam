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
//! * CRAM round-trip — identical mandatory fields (cols 1–11) after decode.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::process::Command;

// ── harness plumbing ──────────────────────────────────────────────────────────

fn samtools() -> Option<String> {
    std::env::var("KIRA_BAM_SAMTOOLS").ok().filter(|s| !s.is_empty())
}

fn kira_bam() -> String {
    env!("CARGO_BIN_EXE_kira-bam").to_string()
}

fn fixtures_dir() -> Option<String> {
    std::env::var("KIRA_BAM_FIXTURES").ok().filter(|s| !s.is_empty())
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
        let s = run_ok("samtools flagstat", Command::new(&st).arg("flagstat").arg(&bam));
        let k = run_ok("kira flagstat", Command::new(kira_bam()).arg("flagstat").arg(&bam));
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
    run_ok("kira sort", Command::new(kira_bam()).args(["sort", "-o"]).arg(&ks).arg(&bam));
    let ss = dir.path().join("s.sorted.bam");
    run_ok("samtools sort", Command::new(&st).args(["sort", "-o"]).arg(&ss).arg(&bam));

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
    run_ok("kira sort", Command::new(kira_bam()).args(["sort", "-o"]).arg(&ks).arg(&bam));
    let ss = dir.path().join("s.bam");
    run_ok("samtools sort", Command::new(&st).args(["sort", "-o"]).arg(&ss).arg(&bam));
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
    run_ok("kira sort -n", Command::new(kira_bam()).args(["sort", "-n", "-o"]).arg(&ks).arg(&bam));
    let ss = dir.path().join("s.bam");
    run_ok("samtools sort -n", Command::new(&st).args(["sort", "-n", "-o"]).arg(&ss).arg(&bam));
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
        Command::new(kira_bam()).args(["merge", "-f", "-o"]).arg(&km).arg(&a).arg(&b),
    );
    let sm = dir.path().join("s.merged.bam");
    run_ok("samtools merge", Command::new(&st).arg("merge").arg("-f").arg(&sm).arg(&a).arg(&b));
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
    run_ok("samtools sort", Command::new(&st).args(["sort", "-o"]).arg(&sorted).arg(&bam));

    // Build the index with kira, then make samtools read it: tests both the
    // kira-written .bai and kira idxstats against samtools.
    run_ok("kira index", Command::new(kira_bam()).arg("index").arg(&sorted));
    let s = run_ok("samtools idxstats", Command::new(&st).arg("idxstats").arg(&sorted));
    let k = run_ok("kira idxstats", Command::new(kira_bam()).arg("idxstats").arg(&sorted));
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
    run_ok("samtools fixmate", Command::new(&st).args(["fixmate", "-m"]).arg(&raw).arg(&fm));
    let prepared = dir.path().join("prepared.bam");
    run_ok("samtools sort", Command::new(&st).args(["sort", "-o"]).arg(&prepared).arg(&fm));

    let s_md = dir.path().join("s.md.bam");
    run_ok("samtools markdup", Command::new(&st).arg("markdup").arg(&prepared).arg(&s_md));
    let k_md = dir.path().join("k.md.bam");
    run_ok("kira markdup", Command::new(kira_bam()).arg("markdup").arg(&prepared).arg(&k_md));

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

    let s = run_ok("samtools calmd", Command::new(&st).arg("calmd").arg(&bam).arg(&refp));
    let k = run_ok("kira calmd", Command::new(kira_bam()).arg("calmd").arg(&bam).arg(&refp));
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
    assert_eq!(original, restored, "CRAM round-trip changed an alignment field");
}

#[test]
fn mpileup_multi_input_columns_match_samtools() {
    let st = require_samtools!();
    let dir = tmpdir();
    let a = build_bam(&st, dir.path(), "ma", MERGE_A_SAM);
    let b = build_bam(&st, dir.path(), "mb", MERGE_B_SAM);

    let k_out = dir.path().join("k.mpileup.tsv");
    run_ok(
        "kira mpileup",
        Command::new(kira_bam())
            .args(["mpileup", "-o"])
            .arg(&k_out)
            .arg(&a)
            .arg(&b),
    );
    let kira = mpileup_depth_projection(&strip_cr(&std::fs::read_to_string(&k_out).unwrap()));

    // -B drops BAQ (kira has none); the projection compares CHROM/POS/REF and
    // every sample's depth only, because kira lacks samtools' ^/$/indel/strand
    // base markup.
    let s_out = run_ok(
        "samtools mpileup",
        Command::new(&st)
            .args(["mpileup", "-B", "-d", "8000", "-q", "0", "-Q", "13"])
            .arg(&a)
            .arg(&b),
    );
    assert_eq!(
        mpileup_depth_projection(&s_out),
        kira,
        "mpileup multi-input projection differs from samtools"
    );
}

// ── comparison helpers ────────────────────────────────────────────────────────

/// CHROM POS REF plus the leading depth field of each following
/// depth/bases/quals triplet — the multi-input invariant both tools must keep.
fn mpileup_depth_projection(out: &str) -> Vec<String> {
    out.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let mut p = vec![f[0].to_string(), f[1].to_string(), f[2].to_string()];
            let mut i = 3;
            while i < f.len() {
                p.push(f[i].to_string());
                i += 3;
            }
            p.join("\t")
        })
        .collect()
}

/// Sorted multiset of non-empty lines (order-insensitive record-set equality).
fn multiset(body: &str) -> Vec<String> {
    let mut v: Vec<String> = body.lines().filter(|l| !l.is_empty()).map(str::to_string).collect();
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
    run_ok("samtools view -f1024", Command::new(samtools).args(["view", "-f", "1024"]).arg(bam))
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
    let header = run_ok("samtools view -H", Command::new(samtools).args(["view", "-H"]).arg(bam));
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
        let tid = if rname == "*" { usize::MAX } else { *order.get(rname).unwrap_or(&usize::MAX) };
        let cur = (tid, pos);
        assert!(cur >= last, "not coordinate-sorted at {line}");
        last = cur;
    }
}
