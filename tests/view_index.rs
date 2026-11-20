use std::process::Command;

#[test]
fn view_uses_an_explicit_index_and_keeps_auxiliary_tags() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    let bam = dir.path().join("input.bam");
    let index = dir.path().join("custom.index");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\nr1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tNM:i:0\n",
    ).unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let output = cli()
        .args(["view", "-u", "-o"])
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    for format in ["-b", "-c"] {
        let output = cli()
            .args(["index", format, "-o"])
            .arg(&index)
            .arg(&bam)
            .output()
            .unwrap();
        assert!(output.status.success(), "{output:?}");
        let output = cli()
            .args(["view", "-X"])
            .arg(&bam)
            .arg(&index)
            .arg("chr1:1-1")
            .output()
            .unwrap();
        assert!(output.status.success(), "{output:?}");
        assert!(String::from_utf8_lossy(&output.stdout).contains("\tNM:i:0"));
    }
    let output = cli()
        .args(["view", "--drop-tags"])
        .arg(&bam)
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    assert!(!String::from_utf8_lossy(&output.stdout).contains("\tNM:"));
}
