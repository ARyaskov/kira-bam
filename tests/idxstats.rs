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

#[test]
fn idxstats_counts_records_when_index_metadata_is_missing() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("reads.sam");
    let bam = dir.path().join("reads.bam");
    std::fs::write(
        &sam,
        concat!(
            "@SQ\tSN:chr1\tLN:10\n",
            "m\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\n",
            "u\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n",
        ),
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
    // BAI magic, one reference, zero bins, zero linear intervals,
    // with both optional record-count fields absent.
    let index = b"BAI\x01\x01\0\0\0\0\0\0\0\0\0\0\0";
    std::fs::write(dir.path().join("reads.bam.bai"), index).unwrap();
    let result = cli().arg("idxstats").arg(&bam).output().unwrap();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(result.stdout, b"chr1\t10\t1\t0\n*\t0\t0\t1\n");
}
