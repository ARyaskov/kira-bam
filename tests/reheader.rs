use std::process::Command;

#[test]
fn reheader_rejects_moved_references_and_allows_positional_renames() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.sam");
    let header = dir.path().join("header.sam");
    let output = dir.path().join("output.bam");
    std::fs::write(
        &input,
        "@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr2\tLN:10\nr\t65\tchr1\t1\t60\t1M\tchr2\t2\t0\tA\tI\n",
    )
    .unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let reheader = || {
        cli()
            .arg("reheader")
            .arg(&header)
            .arg(&input)
            .arg("-o")
            .arg(&output)
            .output()
            .unwrap()
    };

    std::fs::write(&header, "@SQ\tSN:chr2\tLN:10\n@SQ\tSN:chr1\tLN:10\n").unwrap();
    let result = reheader();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("moved from index"));
    assert!(!output.exists());

    std::fs::write(&header, "@SQ\tSN:first\tLN:10\n@SQ\tSN:second\tLN:10\n").unwrap();
    let result = reheader();
    assert!(result.status.success(), "{result:?}");
    let result = cli().arg("view").arg(&output).output().unwrap();
    assert!(result.status.success(), "{result:?}");
    let sam = String::from_utf8(result.stdout).unwrap();
    let fields: Vec<_> = sam.trim_end().split('\t').collect();
    assert_eq!((fields[2], fields[6]), ("first", "second"));
}
