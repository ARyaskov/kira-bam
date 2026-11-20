use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "{output:?}");
    output
}

#[test]
fn view_and_fastq_read_real_sam_and_bgzf_bam_pipes() {
    let sam = b"@HD\tVN:1.6\tSO:queryname\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n";
    let bam = run(&["view", "-b", "-"], sam).stdout;
    for input in [sam.as_slice(), &bam] {
        let count = run(&["view", "-c", "-"], input);
        assert_eq!(String::from_utf8_lossy(&count.stdout).trim(), "1");
        let reads = run(&["fastq", "-"], input);
        assert_eq!(String::from_utf8_lossy(&reads.stdout), "@r1\nA\n+\nI\n");
    }
}
