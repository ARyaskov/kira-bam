use std::collections::HashMap;
use std::process::Command;

#[test]
fn merge_preserves_program_provenance_on_every_record() {
    let dir = tempfile::tempdir().unwrap();
    let first = dir.path().join("first.sam");
    let second = dir.path().join("second.sam");
    let output = dir.path().join("merged.bam");
    for (path, version) in [(&first, 1), (&second, 2)] {
        std::fs::write(
            path,
            format!(
                "@SQ\tSN:chr1\tLN:10\n@PG\tID:aligner\tVN:{version}\n\
                 r{version}a\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tPG:Z:aligner\n\
                 r{version}b\t0\tchr1\t2\t60\t1M\t*\t0\t0\tA\tI\tPG:Z:aligner\n"
            ),
        )
        .unwrap();
    }
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let result = cli()
        .args(["merge", "-o"])
        .arg(&output)
        .args([first, second])
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    let result = cli().args(["view", "-h"]).arg(&output).output().unwrap();
    assert!(result.status.success(), "{result:?}");
    let sam = String::from_utf8(result.stdout).unwrap();
    let versions: HashMap<_, _> = sam
        .lines()
        .filter(|line| line.starts_with("@PG\t"))
        .map(|line| {
            let field = |prefix| {
                line.split('\t')
                    .find_map(|value| value.strip_prefix(prefix))
            };
            (field("ID:").unwrap(), field("VN:"))
        })
        .collect();
    let records: Vec<_> = sam.lines().filter(|line| !line.starts_with('@')).collect();
    assert_eq!(records.len(), 4);
    for record in records {
        let id = record
            .split('\t')
            .find_map(|value| value.strip_prefix("PG:Z:"))
            .unwrap();
        assert_eq!(versions[id], Some(&record[1..2]), "{record}");
    }
}
