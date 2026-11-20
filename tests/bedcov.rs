use std::process::Command;

#[test]
fn bedcov_rejects_invalid_rows_and_reports_zero_for_empty_intervals() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    let bam = dir.path().join("input.bam");
    let bed = dir.path().join("regions.bed");
    std::fs::write(
        &sam,
        "@SQ\tSN:chr1\tLN:10\nr\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\tII\n",
    )
    .unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let result = cli()
        .args(["sort", "-o"])
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    let result = cli().arg("index").arg(&bam).output().unwrap();
    assert!(result.status.success(), "{result:?}");
    let bedcov = || cli().arg("bedcov").arg(&bed).arg(&bam).output().unwrap();

    for row in ["chr1\t0X\t2", "chr1\t0", "chr1\t2\t1"] {
        std::fs::write(&bed, format!("chr1\t0\t2\n{row}\n")).unwrap();
        let result = bedcov();
        assert!(!result.status.success(), "{row}");
        assert!(result.stdout.is_empty(), "{result:?}");
        assert!(String::from_utf8_lossy(&result.stderr).contains("regions.bed: line 2"));
    }
    std::fs::write(&bed, "chr1\t0\t2\nchr1\t1\t1\n").unwrap();
    let result = bedcov();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(result.stdout, b"chr1\t0\t2\t2\nchr1\t1\t1\t0\n");
}
