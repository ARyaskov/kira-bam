use std::process::Command;

#[test]
fn idxstats_rejects_indexes_with_a_different_reference_count() {
    let dir = tempfile::tempdir().unwrap();
    let small = dir.path().join("small.sam");
    let large = dir.path().join("large.sam");
    let small_bam = dir.path().join("small.bam");
    let large_bam = dir.path().join("large.bam");
    std::fs::write(&small, "@SQ\tSN:chr1\tLN:10\n").unwrap();
    std::fs::write(&large, "@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr2\tLN:10\n").unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    for (input, output) in [(&small, &small_bam), (&large, &large_bam)] {
        let result = cli()
            .args(["sort", "-o"])
            .arg(output)
            .arg(input)
            .output()
            .unwrap();
        assert!(result.status.success(), "{result:?}");
    }
    let result = cli().arg("index").arg(&small_bam).output().unwrap();
    assert!(result.status.success(), "{result:?}");
    std::fs::copy(
        dir.path().join("small.bam.bai"),
        dir.path().join("large.bam.bai"),
    )
    .unwrap();
    let result = cli().arg("idxstats").arg(&large_bam).output().unwrap();
    assert!(!result.status.success(), "{result:?}");
    assert!(result.stdout.is_empty());
    assert!(String::from_utf8_lossy(&result.stderr).contains("reference count"));
}
