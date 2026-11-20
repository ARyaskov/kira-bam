use std::process::Command;

#[test]
fn depth_matches_short_quality_flags_and_can_restore_duplicates() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n",
            "lowbq\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\t&\n",
            "lowmq\t0\tchr1\t2\t5\t1M\t*\t0\t0\tT\tI\n",
            "dupe\t1024\tchr1\t2\t60\t1M\t*\t0\t0\tT\tI\n",
        ),
    )
    .unwrap();
    for (flags, expected) in [
        (["-q", "20"], "chr1\t2\t1\n"),
        (["-Q", "20"], "chr1\t1\t1\n"),
        (["-g", "0x400"], "chr1\t1\t1\nchr1\t2\t2\n"),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
            .arg("depth")
            .args(flags)
            .arg(&sam)
            .output()
            .unwrap();
        assert!(result.status.success(), "{result:?}");
        assert_eq!(String::from_utf8_lossy(&result.stdout), expected);
    }
}
