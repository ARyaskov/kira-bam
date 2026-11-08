//! Regression tests for `sort --markdup` memory-spill behaviour:
//! a spilled sort must fail with an actionable error instead of silently
//! emitting sort-only output (no duplicate flags). No samtools needed.

use std::process::Command;

fn kira_bam() -> &'static str {
    env!("CARGO_BIN_EXE_kira-bam")
}

/// A plain SAM file of single-end records big enough that the per-record
/// heap estimate (~670 bytes: 160 fixed + 250 seq + 250 qual + name/cigar)
/// crosses a 1 MiB `-m` budget partway through, forcing a disk spill.
fn spill_sam(records: usize, seq_len: usize) -> String {
    // Cycle the alphabet so SEQ is exactly seq_len bases — CIGAR, SEQ and
    // QUAL must all agree on length or the SAM parser rejects the record.
    let seq: String = (0..seq_len).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect();
    let qual = "I".repeat(seq_len);
    let cigar = format!("{seq_len}M");
    let mut s = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:250000\n");
    for i in 0..records {
        // Spread positions so nothing is dropped as invalid.
        let pos = (i * 17) % 240_000 + 1;
        s.push_str(&format!(
            "r{i}\t0\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n"
        ));
    }
    s
}

#[test]
fn sort_markdup_fails_on_spill_without_output() {
    let dir = tempfile::TempDir::new().unwrap();
    let input = dir.path().join("spill.sam");
    std::fs::write(&input, spill_sam(3000, 250)).unwrap();
    let out = dir.path().join("out.bam");

    let status = Command::new(kira_bam())
        .args(["sort", "--markdup", "-m", "1M", "-o"])
        .arg(&out)
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        !status.status.success(),
        "spilled sort --markdup must fail, not fall back to sort-only"
    );
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(
        stderr.contains("--markdup") && stderr.contains("-m"),
        "stderr must be actionable (name the flag and the fix): {stderr}"
    );
    assert!(
        !out.exists(),
        "a failed sort --markdup must not leave an output file"
    );
}

#[test]
fn sort_markdup_single_chunk_still_succeeds() {
    let dir = tempfile::TempDir::new().unwrap();
    let input = dir.path().join("tiny.sam");
    std::fs::write(&input, spill_sam(10, 10)).unwrap();
    let out = dir.path().join("out.bam");

    let status = Command::new(kira_bam())
        .args(["sort", "--markdup", "-o"])
        .arg(&out)
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        status.status.success(),
        "single-chunk sort --markdup must keep working: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    assert!(out.exists());
}
