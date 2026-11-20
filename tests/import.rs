use std::process::Command;

#[test]
fn paired_import_validates_both_mates() {
    let dir = tempfile::tempdir().unwrap();
    let r1 = dir.path().join("r1.fq");
    let r2 = dir.path().join("r2.fq");
    let bam = dir.path().join("out.bam");
    let record1 = "@pair/1\nA\n+\nI\n";
    let record2 = "@pair/2\nT\n+\nI\n";
    let run = |left, right| {
        std::fs::write(&r1, left).unwrap();
        std::fs::write(&r2, right).unwrap();
        Command::new(env!("CARGO_BIN_EXE_kira-bam"))
            .args(["import", "-o"])
            .arg(&bam)
            .arg(&r1)
            .arg(&r2)
            .output()
            .unwrap()
    };

    for (left, right, message) in [
        (record1, "", "different record counts"),
        ("", record2, "different record counts"),
        (record1, "@shifted/2\nT\n+\nI\n", "mate names differ"),
        (record1, "@pair/2\nT\n+\n", "truncated FASTQ record"),
        (record1, "@pair/2\nT\n+\n\n", "lengths differ"),
    ] {
        let output = run(left, right);
        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(message),
            "{output:?}"
        );
    }
    let output = run(record1, record2);
    assert!(output.status.success(), "{output:?}");
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(["view", "-c"])
        .arg(&bam)
        .output()
        .unwrap();
    assert_eq!(String::from_utf8_lossy(&output.stdout).trim(), "2");
}
