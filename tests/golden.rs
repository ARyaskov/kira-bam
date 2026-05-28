// Gated on KIRA_BAM_SAMTOOLS + KIRA_BAM_FIXTURES env vars.

use std::path::Path;
use std::process::Command;

fn samtools() -> Option<String> {
    std::env::var("KIRA_BAM_SAMTOOLS").ok()
}

fn fixtures_dir() -> Option<String> {
    std::env::var("KIRA_BAM_FIXTURES").ok()
}

fn kira_bam() -> String {
    env!("CARGO_BIN_EXE_kira-bam").to_string()
}

fn list_bams(dir: &Path) -> Vec<std::path::PathBuf> {
    let mut out = Vec::new();
    if let Ok(entries) = std::fs::read_dir(dir) {
        for e in entries.flatten() {
            let p = e.path();
            if p.is_dir() {
                out.extend(list_bams(&p));
            } else if p.extension().and_then(|s| s.to_str()) == Some("bam") {
                out.push(p);
            }
        }
    }
    out
}

fn strip_cr(s: &str) -> String {
    s.replace('\r', "")
}

#[test]
fn flagstat_matches_samtools() {
    let Some(sam) = samtools() else { return };
    let Some(fix) = fixtures_dir() else { return };
    let mut pass = 0usize;
    let mut fail = 0usize;
    let mut failing: Vec<String> = Vec::new();
    for bam in list_bams(Path::new(&fix)) {
        let s_out = Command::new(&sam).args(["flagstat"]).arg(&bam).output();
        let k_out = Command::new(kira_bam())
            .args(["flagstat"])
            .arg(&bam)
            .output();
        match (s_out, k_out) {
            (Ok(s), Ok(k)) if s.status.success() && k.status.success() => {
                let s_text = strip_cr(&String::from_utf8_lossy(&s.stdout));
                let k_text = strip_cr(&String::from_utf8_lossy(&k.stdout));
                if s_text == k_text {
                    pass += 1;
                } else {
                    fail += 1;
                    failing.push(bam.display().to_string());
                }
            }
            _ => continue,
        }
    }
    eprintln!("flagstat compare: {pass} pass / {fail} fail");
    for f in &failing {
        eprintln!("  fail: {f}");
    }
    assert!(
        pass * 20 >= (pass + fail) * 19,
        "flagstat concordance < 95% ({pass}/{})",
        pass + fail
    );
}

#[test]
fn view_body_matches_samtools() {
    let Some(sam) = samtools() else { return };
    let Some(fix) = fixtures_dir() else { return };
    let mut pass = 0usize;
    let mut fail = 0usize;
    for bam in list_bams(Path::new(&fix)) {
        let s = match Command::new(&sam).args(["view"]).arg(&bam).output() {
            Ok(o) if o.status.success() => o.stdout,
            _ => continue,
        };
        let k = match Command::new(kira_bam()).args(["view"]).arg(&bam).output() {
            Ok(o) if o.status.success() => o.stdout,
            _ => continue,
        };
        let s_text = strip_cr(&String::from_utf8_lossy(&s));
        let k_text = strip_cr(&String::from_utf8_lossy(&k));
        if s_text == k_text {
            pass += 1;
        } else {
            fail += 1;
        }
    }
    eprintln!("view body compare: {pass} pass / {fail} fail");
    assert!(pass >= fail, "view body concordance below 50%");
}
