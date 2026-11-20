use std::process::Command;

#[test]
fn fqidx_preserves_crlf_offsets_and_rejects_an_incomplete_record() {
    let dir = tempfile::tempdir().unwrap();
    let fastq = dir.path().join("reads.fq");
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    std::fs::write(&fastq, "@read/1 comment\r\nACT\r\n+\r\nIJK\r\n").unwrap();
    let output = cli()
        .arg("fqidx")
        .arg(&fastq)
        .arg("read/1")
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    assert_eq!(output.stdout, b"@read/1\nACT\n+\nIJK\n");

    let truncated = dir.path().join("truncated.fq");
    std::fs::write(&truncated, "@read\nA\n+\n").unwrap();
    let output = cli().arg("fqidx").arg(&truncated).output().unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("truncated FASTQ record"));
    assert!(!dir.path().join("truncated.fq.fqi").exists());
}
