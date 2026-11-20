use std::process::Command;

#[test]
fn fastq_routes_paired_and_unpaired_read2_to_its_own_stream() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    let reads = dir.path().join("r2.fq");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:queryname\n",
            "pair\t77\t*\t0\t0\t*\t*\t0\t0\tA\tI\n",
            "pair\t141\t*\t0\t0\t*\t*\t0\t0\tT\tI\n",
            "solo\t141\t*\t0\t0\t*\t*\t0\t0\tC\tI\n",
        ),
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(["fastq", "-2"])
        .arg(&reads)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    assert_eq!(String::from_utf8_lossy(&output.stdout), "@pair\nA\n+\nI\n");
    assert_eq!(
        std::fs::read_to_string(reads).unwrap(),
        "@pair\nT\n+\nI\n@solo\nC\n+\nI\n"
    );
}
