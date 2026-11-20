use std::process::Command;

#[test]
fn collate_enables_fixmate_and_coordinate_sort_clears_query_grouping() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.sam");
    let grouped = dir.path().join("grouped.bam");
    let fixed = dir.path().join("fixed.bam");
    let sorted = dir.path().join("sorted.bam");
    std::fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\tSS:coordinate:minhash\n",
            "@SQ\tSN:chr1\tLN:100\n",
            "pair\t65\tchr1\t1\t60\t1M\t=\t10\t0\tA\tI\n",
            "other\t0\tchr1\t5\t60\t1M\t*\t0\t0\tG\tI\n",
            "pair\t145\tchr1\t10\t60\t1M\t=\t1\t0\tT\tI\n",
        ),
    )
    .unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let output = cli()
        .arg("collate")
        .arg(&input)
        .arg(&grouped)
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let output = cli()
        .args(["fixmate", "-m"])
        .arg(&grouped)
        .arg(&fixed)
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let output = cli().args(["view", "-h"]).arg(&fixed).output().unwrap();
    let sam = std::str::from_utf8(&output.stdout).unwrap();
    assert!(sam.contains("SO:unsorted\tGO:query"));
    assert!(!sam.contains("SS:coordinate:minhash"));
    assert_eq!(
        sam.lines()
            .filter(|line| line.starts_with("pair\t") && line.contains("MC:Z:1M"))
            .count(),
        2
    );

    let output = cli()
        .args(["sort", "-o"])
        .arg(&sorted)
        .arg(&grouped)
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let output = cli()
        .arg("fixmate")
        .arg(&sorted)
        .arg(&fixed)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("needs the two mates"));
}
