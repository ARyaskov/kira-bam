use std::process::Command;

#[test]
fn fqidx_rebuilds_stale_offsets_and_rejects_malformed_fields() {
    let dir = tempfile::tempdir().unwrap();
    let fastq = dir.path().join("reads.fq");
    let index = dir.path().join("reads.fq.fqi");
    let cli = || {
        Command::new(env!("CARGO_BIN_EXE_kira-bam"))
            .arg("fqidx")
            .arg(&fastq)
            .arg("a")
            .output()
            .unwrap()
    };
    std::fs::write(&fastq, "@a\nAAA\n+\nIII\n@b\nCCC\n+\nJJJ\n").unwrap();
    let result = cli();
    assert!(result.status.success(), "{result:?}");

    std::fs::write(&fastq, "@b\nCCC\n+\nJJJ\n@a\nAAA\n+\nIII\n").unwrap();
    std::fs::OpenOptions::new()
        .write(true)
        .open(&index)
        .unwrap()
        .set_modified(std::time::UNIX_EPOCH)
        .unwrap();
    let result = cli();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(result.stdout, b"@a\nAAA\n+\nIII\n");

    for (entry, message) in [
        ("a\t3\tbad\t8\n", "invalid sequence offset"),
        ("a\t3\t4\t9999\n", "inconsistent"),
        ("a\t3\n", "expected four"),
    ] {
        std::fs::write(&index, entry).unwrap();
        let result = cli();
        assert!(!result.status.success());
        assert!(result.stdout.is_empty());
        let error = String::from_utf8_lossy(&result.stderr);
        assert!(
            error.contains("line 1") && error.contains(message),
            "{error}"
        );
    }
}

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

#[test]
fn fqidx_rejects_output_aliases_before_creating_the_index() {
    let dir = tempfile::tempdir().unwrap();
    let fastq = dir.path().join("input.fq");
    let index = dir.path().join("input.fq.fqi");
    let reads = b"@read\nAC\n+\nII\n";
    std::fs::write(&fastq, reads).unwrap();
    for destination in [&fastq, &index] {
        let result = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
            .arg("fqidx")
            .arg(&fastq)
            .arg("read")
            .arg("-o")
            .arg(destination)
            .output()
            .unwrap();
        assert!(!result.status.success(), "{result:?}");
        assert_eq!(std::fs::read(&fastq).unwrap(), reads);
        assert!(!index.exists());
    }
}
