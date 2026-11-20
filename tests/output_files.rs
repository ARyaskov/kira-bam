use std::process::Command;

#[test]
fn depth_keeps_previous_report_after_a_late_input_error() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("bad-order.sam");
    let report = dir.path().join("depth.tsv");
    std::fs::write(&report, "previous report\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n",
            "r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\n",
            "r2\t0\tchr1\t6\t60\t1M\t*\t0\t0\tT\tI\n",
            "r3\t0\tchr1\t2\t60\t1M\t*\t0\t0\tT\tI\n",
        ),
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(["depth", "-o"])
        .arg(&report)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("not coordinate-sorted"));
    assert_eq!(
        std::fs::read_to_string(report).unwrap(),
        "previous report\n"
    );
}

#[test]
fn reheader_rejects_a_hard_link_to_its_input() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("source.sam");
    let alias = dir.path().join("alias.bam");
    let header = dir.path().join("header.sam");
    let data = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    std::fs::write(&sam, data).unwrap();
    std::fs::write(&header, data).unwrap();
    std::fs::hard_link(&sam, &alias).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(["reheader", "-o"])
        .arg(&alias)
        .arg(&header)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("refers to input"));
    assert_eq!(std::fs::read_to_string(sam).unwrap(), data);
}

#[test]
fn split_collision_preserves_a_previous_result() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("source.sam");
    let bam = dir.path().join("result.bam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\n",
            "r1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\tRG:Z:a\n",
            "r2\t4\t*\t0\t0\t*\t*\t0\t0\tT\tI\tRG:Z:b\n",
        ),
    )
    .unwrap();
    std::fs::write(&bam, b"previous output").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(["split", "-f"])
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("multiple tag values"));
    assert_eq!(std::fs::read(bam).unwrap(), b"previous output");
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
}
