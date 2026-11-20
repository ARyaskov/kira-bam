use std::process::Command;

#[test]
fn coverage_clips_bases_and_denominator_to_the_requested_interval() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n",
            "left\t0\tchr1\t8\t60\t4M\t*\t0\t0\tAAAA\tIIII\n",
            "dup\t1024\tchr1\t10\t60\t4M\t*\t0\t0\tAAAA\tIIII\n",
            "right\t0\tchr1\t13\t40\t4M\t*\t0\t0\tAAAA\t5555\n",
        ),
    )
    .unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(["coverage", "-H", "-r", "chr1:10-13"])
        .arg(&sam)
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(
        String::from_utf8_lossy(&result.stdout),
        "chr1\t10\t13\t2\t3\t75.0000\t0.7500\t33.3333\t50.0000\n",
    );
}

#[test]
fn indexed_coverage_matches_streaming_over_overlapping_and_disjoint_regions() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    let bam = dir.path().join("input.bam");
    std::fs::write(
        &sam,
        concat!(
            "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:100\n",
            "long\t0\tchr1\t1\t60\t12M\t*\t0\t0\tAAAAAAAAAAAA\tIIIIIIIIIIII\n",
            "other\t0\tchr2\t2\t30\t2M\t*\t0\t0\tAA\tII\n",
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
    let result = cli().arg("index").arg(&bam).output().unwrap();
    assert!(result.status.success(), "{result:?}");
    let coverage = |path: &std::path::Path| {
        let result = cli()
            .args([
                "coverage",
                "-r",
                "chr1:9-12",
                "-r",
                "chr1:3-5",
                "-r",
                "chr1:2-4",
                "-r",
                "chr2:1-3",
            ])
            .arg(path)
            .output()
            .unwrap();
        assert!(result.status.success(), "{result:?}");
        result.stdout
    };
    assert_eq!(coverage(&bam), coverage(&sam));
}
