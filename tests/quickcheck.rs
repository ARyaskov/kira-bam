use kira_bam::cli::QuickcheckArgs;
use kira_bam::quickcheck::{CheckFailures, run};
use std::process::Command;

#[test]
fn quickcheck_returns_file_errors_to_library_callers() {
    let dir = tempfile::tempdir().unwrap();
    let missing = dir.path().join("missing.bam");
    let err = run(QuickcheckArgs {
        inputs: vec![missing.clone()],
        verbose: false,
        unrecognized_ok: false,
    })
    .unwrap_err();
    let failures = err.downcast_ref::<CheckFailures>().unwrap();
    assert_eq!(failures.files[0].0, missing);

    let output = Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(["quickcheck", "-v"])
        .arg(&missing)
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&output.stderr).contains("missing.bam"));
}
