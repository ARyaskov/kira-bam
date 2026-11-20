use std::process::Command;

#[test]
fn split_distinguishes_integer_text_and_missing_tag_values() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.sam");
    let mut sam = String::from("@SQ\tSN:chr1\tLN:10\n");
    for (name, tag) in [
        ("zero", "\tXX:i:0"),
        ("one", "\tXX:i:1"),
        ("text_zero", "\tXX:Z:0"),
        ("missing", ""),
        ("literal", "\tXX:Z:unknown"),
    ] {
        sam.push_str(&format!("{name}\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI{tag}\n"));
    }
    std::fs::write(&input, sam).unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let result = cli()
        .current_dir(dir.path())
        .args(["split", "-d", "XX"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    for (suffix, name) in [
        ("%i0", "zero"),
        ("%i1", "one"),
        ("0", "text_zero"),
        ("unknown", "missing"),
        ("%75nknown", "literal"),
    ] {
        let result = cli()
            .arg("view")
            .arg(dir.path().join(format!("input_{suffix}.bam")))
            .output()
            .unwrap();
        assert!(result.status.success(), "{result:?}");
        let sam = String::from_utf8(result.stdout).unwrap();
        let names: Vec<_> = sam
            .lines()
            .map(|line| line.split('\t').next().unwrap())
            .collect();
        assert_eq!(names, [name]);
    }
}
