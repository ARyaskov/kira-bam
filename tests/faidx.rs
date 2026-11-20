use std::process::Command;

#[test]
fn faidx_guards_source_paths_and_validates_index_geometry() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("reference.fa");
    let index = dir.path().join("reference.fa.fai");
    let link = dir.path().join("alias.fa");
    let original = b">chr1\r\nACGT\r\nTT\r\n";
    std::fs::write(&fasta, original).unwrap();
    std::fs::hard_link(&fasta, &link).unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));

    for destination in [&fasta, &index, &link] {
        let result = cli()
            .arg("faidx")
            .arg(&fasta)
            .arg("chr1:2-5")
            .arg("-o")
            .arg(destination)
            .output()
            .unwrap();
        assert!(!result.status.success(), "{result:?}");
        assert_eq!(std::fs::read(&fasta).unwrap(), original);
        assert!(!index.exists());
    }
    let query = || {
        cli()
            .arg("faidx")
            .arg(&fasta)
            .arg("chr1:2-5")
            .output()
            .unwrap()
    };
    let result = query();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(result.stdout, b">chr1:2-5\nCGTT\n");

    std::fs::write(&index, "chr1\t6\tinvalid\t4\t6\n").unwrap();
    let result = query();
    assert!(!result.status.success(), "{result:?}");
    assert!(String::from_utf8_lossy(&result.stderr).contains("line 1"));

    let invalid = dir.path().join("invalid.fa");
    std::fs::write(&invalid, b">chr1\nAAAA\nAA\nAAAA\n").unwrap();
    let result = cli().arg("faidx").arg(&invalid).output().unwrap();
    assert!(!result.status.success(), "{result:?}");
    assert!(String::from_utf8_lossy(&result.stderr).contains("line 4"));
}

#[test]
fn faidx_index_build_reports_write_failure_but_queries_can_use_in_memory_offsets() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("reference.fa");
    let index = dir.path().join("reference.fa.fai");
    std::fs::write(&fasta, b">chr1\nACGT\n").unwrap();
    std::fs::create_dir(&index).unwrap();
    // Force a rebuild even on filesystems whose timestamps have coarse precision.
    std::fs::OpenOptions::new()
        .write(true)
        .open(&fasta)
        .unwrap()
        .set_modified(std::time::SystemTime::now() + std::time::Duration::from_secs(60))
        .unwrap();
    let cli = || {
        let mut command = Command::new(env!("CARGO_BIN_EXE_kira-bam"));
        command.arg("faidx").arg(&fasta);
        command
    };
    let result = cli().output().unwrap();
    assert!(!result.status.success(), "{result:?}");
    assert!(String::from_utf8_lossy(&result.stderr).contains("write"));
    let result = cli().arg("chr1:1-2").output().unwrap();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(result.stdout, b">chr1:1-2\nAC\n");
    assert!(index.is_dir());
}
