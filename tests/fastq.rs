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
        std::fs::read_to_string(&reads).unwrap(),
        "@pair\nT\n+\nI\n@solo\nC\n+\nI\n"
    );
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .arg("fastq")
        .arg("-1")
        .arg(&reads)
        .arg("-2")
        .arg(&reads)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    assert!(output.stdout.is_empty());
    assert_eq!(
        std::fs::read_to_string(reads).unwrap(),
        "@pair\nA\n+\nI\n@pair\nT\n+\nI\n@solo\nC\n+\nI\n",
    );
}

#[test]
fn collate_and_fastq_can_stream_a_name_group_larger_than_the_budget() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    let bam = dir.path().join("output.bam");
    let seq = "A".repeat(4096);
    let qual = "I".repeat(4096);
    let count = 400;
    let mut data = String::from("@SQ\tSN:chr1\tLN:10000\n");
    for i in 0..count {
        let flags = if i % 2 == 0 { 77 } else { 141 };
        data.push_str(&format!(
            "repeated\t{flags}\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\n"
        ));
    }
    std::fs::write(&sam, data).unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let result = cli()
        .args(["collate", "-m", "1M", "--buckets", "1"])
        .arg(&sam)
        .arg(&bam)
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    let result = cli().args(["view", "-c"]).arg(&bam).output().unwrap();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(
        String::from_utf8(result.stdout).unwrap().trim(),
        count.to_string()
    );

    let result = cli()
        .args(["fastq", "--collate", "-m", "1M", "--buckets", "1", "-N"])
        .arg(&sam)
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(
        String::from_utf8(result.stdout).unwrap(),
        format!("@repeated/1\n{seq}\n+\n{qual}\n@repeated/2\n{seq}\n+\n{qual}\n")
    );
}
