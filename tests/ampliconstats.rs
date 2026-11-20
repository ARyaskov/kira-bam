use std::process::Command;

#[test]
fn ampliconstats_distinguishes_empty_intervals_from_missing_indexes() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    let bam = dir.path().join("input.bam");
    let bed = dir.path().join("targets.bed");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\nr1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\n",
    )
    .unwrap();
    std::fs::write(&bed, "# targets\nchr1\t0\t1\thit\nchr1\t9\t10\tempty\n").unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    assert!(
        cli()
            .args(["view", "-b", "-o"])
            .arg(&bam)
            .arg(&sam)
            .status()
            .unwrap()
            .success()
    );

    let query = || {
        cli()
            .arg("ampliconstats")
            .arg(&bed)
            .arg(&bam)
            .output()
            .unwrap()
    };
    let unindexed = query();
    assert!(!unindexed.status.success());
    assert!(unindexed.stdout.is_empty());
    assert!(String::from_utf8_lossy(&unindexed.stderr).contains("input.bam"));

    assert!(cli().arg("index").arg(&bam).status().unwrap().success());
    let indexed = query();
    assert!(indexed.status.success(), "{:?}", indexed);
    assert_eq!(
        std::str::from_utf8(&indexed.stdout).unwrap(),
        "# amplicon\tchrom\tstart\tend\tlength\treads...\nhit\tchr1\t1\t1\t1\t1\nempty\tchr1\t10\t10\t1\t0\n"
    );
}
