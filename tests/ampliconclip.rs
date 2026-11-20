use std::process::Command;

#[test]
fn ampliconclip_clears_stale_alignment_tags_and_declares_unsorted_output() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    let bam = dir.path().join("output.bam");
    let bed = dir.path().join("primers.bed");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n",
            "full\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\tII\tMD:Z:2\tNM:i:0\tMC:Z:2M\n",
            "part\t0\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII\tMD:Z:4\tNM:i:0\n",
        ),
    )
    .unwrap();
    std::fs::write(&bed, "chr1\t0\t2\n").unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let result = cli()
        .args(["ampliconclip", "-b"])
        .arg(&bed)
        .arg("-o")
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    let result = cli().args(["view", "-h"]).arg(&bam).output().unwrap();
    assert!(result.status.success(), "{result:?}");
    let output = String::from_utf8(result.stdout).unwrap();
    assert!(output.lines().next().unwrap().contains("SO:unsorted"));
    let full = output
        .lines()
        .find(|line| line.starts_with("full\t"))
        .unwrap();
    assert_eq!(full.split('\t').nth(1), Some("4"));
    for tag in ["\tMD:", "\tNM:", "\tMC:"] {
        assert!(!output.contains(tag), "{output}");
    }
    assert!(output.contains("part\t0\tchr1\t3\t60\t2S2M"));
}
