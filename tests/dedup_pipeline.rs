use std::collections::BTreeSet;
use std::process::{Command, Output};

fn checked(command: &mut Command) -> Output {
    let output = command.output().unwrap();
    assert!(output.status.success(), "{command:?}: {output:?}");
    output
}

fn duplicate_reads(sam: &[u8]) -> BTreeSet<(String, u16)> {
    std::str::from_utf8(sam)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('@'))
        .filter_map(|line| {
            let (name, rest) = line.split_once('\t').unwrap();
            let flags: u16 = rest.split('\t').next().unwrap().parse().unwrap();
            (flags & 0x400 != 0).then(|| (name.to_owned(), flags & 0xc0))
        })
        .collect()
}

#[test]
fn collate_fixmate_sort_markdup_preserve_both_mates_of_each_template() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("reads.sam");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n",
            "winner\t99\tchr1\t1\t60\t4M\t=\t12\t15\tACGT\tIIII\n",
            "loser\t99\tchr1\t1\t60\t4M\t=\t12\t15\tACGT\t5555\n",
            "winner\t147\tchr1\t12\t60\t4M\t=\t1\t-15\tACGT\tIIII\n",
            "loser\t147\tchr1\t12\t60\t4M\t=\t1\t-15\tACGT\t5555\n",
        ),
    )
    .unwrap();
    let cli = || Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    let grouped = dir.path().join("grouped.bam");
    let fixed = dir.path().join("fixed.bam");
    let sorted = dir.path().join("sorted.bam");
    let marked = dir.path().join("marked.bam");
    checked(cli().arg("collate").arg(&sam).arg(&grouped));
    checked(cli().args(["fixmate", "-m"]).arg(&grouped).arg(&fixed));
    checked(cli().args(["sort", "-o"]).arg(&sorted).arg(&fixed));
    checked(cli().arg("markdup").arg(&sorted).arg(&marked));
    let decoded = checked(cli().arg("view").arg(&marked)).stdout;
    assert_eq!(std::str::from_utf8(&decoded).unwrap().lines().count(), 4);
    let actual = duplicate_reads(&decoded);
    assert_eq!(
        actual,
        BTreeSet::from([(String::from("loser"), 0x40), (String::from("loser"), 0x80)]),
    );

    if let Some(samtools) = std::env::var_os("KIRA_BAM_SAMTOOLS") {
        let cli = || Command::new(&samtools);
        let raw = dir.path().join("oracle.raw.bam");
        let prefix = dir.path().join("oracle.tmp");
        let grouped = dir.path().join("oracle.grouped.bam");
        let fixed = dir.path().join("oracle.fixed.bam");
        let sorted = dir.path().join("oracle.sorted.bam");
        let marked = dir.path().join("oracle.marked.bam");
        checked(cli().args(["view", "-b", "-o"]).arg(&raw).arg(&sam));
        checked(
            cli()
                .args(["collate", "-o"])
                .arg(&grouped)
                .arg(&raw)
                .arg(prefix),
        );
        checked(cli().args(["fixmate", "-m"]).arg(&grouped).arg(&fixed));
        checked(cli().args(["sort", "-o"]).arg(&sorted).arg(&fixed));
        checked(cli().arg("markdup").arg(&sorted).arg(&marked));
        assert_eq!(
            actual,
            duplicate_reads(&checked(cli().arg("view").arg(&marked)).stdout)
        );
    }
}
