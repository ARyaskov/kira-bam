//! One differential test of the mpileup markup against a real samtools.
//!
//! Skipped — with a distinguishable `OPENCODE_ORACLE` receipt — unless a
//! samtools is reachable: `KIRA_BAM_SAMTOOLS` (the golden tests' env var) or a
//! `wsl -e samtools` lead. A tiny SAM with insertions, deletions, a reference
//! skip, reverse-strand reads, a MAPQ-255 read, an all-filtered `-Q` read and
//! a `=`-nibble read runs through both tools with no reference; the rows must
//! match byte for byte. Expected values come from the samtools contract, not
//! from kira's own output.

use std::path::Path;
use std::process::Command;

const SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
f1\t0\tchr1\t11\t30\t4M1I3M1D2M3N2M\t*\t0\t0\tAAAACCCCTTAA\tIIIIIIIIIIIA
r1\t16\tchr1\t11\t20\t5M\t*\t0\t0\tGGGGG\tIIIII
r2\t16\tchr1\t11\t30\t2M2D2M1I2M\t*\t0\t0\tAACCCAA\tIIIIIII
sec\t260\tchr1\t11\t30\t5M\t*\t0\t0\tAAAAA\tIIIII
sup\t2048\tchr1\t40\t30\t5M\t*\t0\t0\tAAAAA\tIIIII
unm\t4\tchr1\t50\t30\t5M\t*\t0\t0\tAAAAA\tIIIII
lowq\t0\tchr1\t60\t30\t5M\t*\t0\t0\tAAAAA\t!!!!!
eq\t0\tchr1\t70\t30\t5=\t*\t0\t0\t=====\tIIIII
mq255\t0\tchr1\t80\t255\t5M\t*\t0\t0\tAAAAA\tIIIII
";

fn strip_cr(s: &str) -> String {
    s.replace('\r', "")
}

fn run_ok(label: &str, mut cmd: Command) -> String {
    let out = cmd
        .output()
        .unwrap_or_else(|e| panic!("{label}: spawn failed: {e}"));
    assert!(
        out.status.success(),
        "{label}: exit {:?}\nstderr: {}",
        out.status.code(),
        String::from_utf8_lossy(&out.stderr)
    );
    strip_cr(&String::from_utf8_lossy(&out.stdout))
}

fn first_line(s: &str) -> String {
    s.lines().next().unwrap_or("").trim().to_string()
}

/// Escape a string for embedding in a JSON string literal.
fn json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

/// A resolved samtools oracle: command factory, path converter, resolved
/// binary path and observed `--version` line.
type RunFactory = Box<dyn Fn(&[&str]) -> Command>;
type PathConverter = Box<dyn Fn(&Path) -> String>;

struct Oracle {
    run: RunFactory,
    to_os: PathConverter,
    path: String,
    version: String,
}

/// Locate a usable samtools: the golden tests' env var first, then the WSL
/// lead. The version is observed from a real `--version` run, never
/// hardcoded. `None` = genuine absence (the caller emits a SKIPPED receipt).
fn probe_oracle() -> Option<Oracle> {
    if let Ok(s) = std::env::var("KIRA_BAM_SAMTOOLS")
        && !s.is_empty()
    {
        let out = Command::new(&s).arg("--version").output().ok()?;
        if !out.status.success() {
            return None;
        }
        let version = first_line(&String::from_utf8_lossy(&out.stdout));
        let path = s.clone();
        return Some(Oracle {
            run: Box::new(move |argv: &[&str]| {
                let mut c = Command::new(&s);
                c.args(argv);
                c
            }),
            to_os: Box::new(|p| p.to_string_lossy().into_owned()),
            path,
            version,
        });
    }
    // WSL lead: `wsl -e samtools --version` doubles as the presence probe.
    let mut ver = Command::new("wsl");
    ver.args(["-e", "samtools", "--version"]);
    let out = ver.output().ok()?;
    if !out.status.success() {
        return None;
    }
    let version = first_line(&String::from_utf8_lossy(&out.stdout));
    // Resolve the binary's absolute path inside WSL for the receipt.
    let mut which = Command::new("wsl");
    which.args(["-e", "which", "samtools"]);
    let path = which
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| first_line(&String::from_utf8_lossy(&o.stdout)))
        .filter(|p| !p.is_empty())
        .unwrap_or_else(|| "wsl:samtools".to_string());
    Some(Oracle {
        run: Box::new(|argv: &[&str]| {
            let mut c = Command::new("wsl");
            c.arg("-e").arg("samtools").args(argv);
            c
        }),
        to_os: Box::new(wsl_path),
        path,
        version,
    })
}

/// WSL-visible form of a Windows path (`C:\x\y` → `/mnt/c/x/y`).
fn wsl_path(p: &Path) -> String {
    let s = p.to_string_lossy().replace('\\', "/");
    let b = s.as_bytes();
    if b.len() >= 2 && b[1] == b':' {
        format!("/mnt/{}{}", (b[0] as char).to_ascii_lowercase(), &s[2..])
    } else {
        s
    }
}

#[test]
fn mpileup_markup_matches_samtools() {
    let Some(oracle) = probe_oracle() else {
        println!(
            "OPENCODE_ORACLE {{\"test\":\"mpileup_markup_vs_samtools\",\"status\":\"SKIPPED\",\"reason\":\"no samtools on KIRA_BAM_SAMTOOLS/PATH/WSL\"}}"
        );
        return;
    };

    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, SAM).unwrap();
    let bam = dir.path().join("in.bam");

    let mut build = (oracle.run)(&["view", "-b", "-o"]);
    build.arg((oracle.to_os)(&bam)).arg((oracle.to_os)(&sam));
    run_ok("build_bam", build);

    let kira_out = dir.path().join("kira.tsv");
    let mut k = Command::new(env!("CARGO_BIN_EXE_kira-bam"));
    k.args(["mpileup", "-o"]).arg(&kira_out).arg(&bam);
    run_ok("kira mpileup", k);
    let kira = strip_cr(&std::fs::read_to_string(&kira_out).unwrap());

    // samtools writes the pileup to stdout; its stderr carries progress
    // lines, which run_ok discards.
    let mut s = (oracle.run)(&["mpileup"]);
    s.arg((oracle.to_os)(&bam));
    let sam = run_ok("samtools mpileup", s);

    assert_eq!(kira, sam, "kira mpileup rows diverge from samtools");

    // Receipt only after the real subprocess ran and the comparison passed.
    println!(
        "OPENCODE_ORACLE {{\"test\":\"mpileup_markup_vs_samtools\",\"status\":\"EXECUTED\",\"path\":\"{}\",\"version\":\"{}\",\"oracle_exit_code\":0,\"comparison_passed\":true}}",
        json(&oracle.path),
        json(&oracle.version),
    );
}

#[test]
fn mpileup_adjusts_matching_and_disagreeing_overlap_qualities() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("pair.sam");
    std::fs::write(
        &sam,
        concat!(
            "@SQ\tSN:chr1\tLN:10\n",
            "pair\t99\tchr1\t1\t60\t3M\t=\t2\t4\tAAC\t555\n",
            "pair\t147\tchr1\t2\t60\t3M\t=\t1\t-4\tAGC\t???\n",
        ),
    )
    .unwrap();
    let pileup = |extra: &[&str]| {
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_kira-bam"))
            .arg("mpileup")
            .args(extra)
            .arg(&sam)
            .output()
            .unwrap();
        assert!(result.status.success(), "{result:?}");
        String::from_utf8(result.stdout).unwrap()
    };
    let adjusted = pileup(&[]);
    let rows: Vec<_> = adjusted
        .lines()
        .map(|row| row.split('\t').collect::<Vec<_>>())
        .collect();
    assert_eq!((rows[1][3], rows[1][5]), ("1", "S"));
    assert_eq!((rows[2][3], rows[2][4], rows[2][5]), ("1", "g", "9"));

    let unadjusted = pileup(&["-x"]);
    assert_eq!(
        unadjusted.lines().nth(1).unwrap().split('\t').nth(3),
        Some("2")
    );
    if let Some(samtools) = std::env::var_os("KIRA_BAM_SAMTOOLS") {
        let result = std::process::Command::new(samtools)
            .arg("mpileup")
            .arg(&sam)
            .output()
            .unwrap();
        assert!(result.status.success(), "{result:?}");
        assert_eq!(result.stdout, adjusted.as_bytes());
    }
}

#[test]
fn mpileup_baq_mode_matches_samtools() {
    let Some(samtools) = std::env::var_os("KIRA_BAM_SAMTOOLS") else {
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("ref.fa");
    let sam = dir.path().join("input.sam");
    std::fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();
    std::fs::write(
        &sam,
        concat!(
            "@SQ\tSN:chr1\tLN:16\n",
            "r\t0\tchr1\t2\t60\t2M1I3M\t*\t0\t0\tCGATAA\tIIIIII\n",
        ),
    )
    .unwrap();
    let expected = std::process::Command::new(&samtools)
        .args(["mpileup", "-f"])
        .arg(&fasta)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(expected.status.success(), "{expected:?}");
    assert!(!expected.stdout.is_empty());
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_kira-bam"))
        .args(["mpileup", "--baq", "-f"])
        .arg(&fasta)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    assert_eq!(result.stdout, expected.stdout);
}
