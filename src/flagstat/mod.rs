use anyhow::{Context, Result};

use crate::cli::FlagstatArgs;
use crate::io::{BamReader, RecordSummary, install_thread_pool};

const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAP: u16 = 0x4;
const FLAG_MUNMAP: u16 = 0x8;
const FLAG_READ1: u16 = 0x40;
const FLAG_READ2: u16 = 0x80;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QCFAIL: u16 = 0x200;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

#[derive(Default, Clone, Copy, Debug)]
pub struct Counts {
    pub total: u64,
    pub primary: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub duplicates: u64,
    pub primary_duplicates: u64,
    pub mapped: u64,
    pub primary_mapped: u64,
    pub paired: u64,
    pub read1: u64,
    pub read2: u64,
    pub properly_paired: u64,
    pub with_itself_and_mate_mapped: u64,
    pub singletons: u64,
    pub mate_diff_chr: u64,
    pub mate_diff_chr_mq5: u64,
}

#[derive(Default, Clone, Copy, Debug)]
pub struct FlagstatTotals {
    pub pass: Counts,
    pub fail: Counts,
}

impl FlagstatTotals {
    pub fn update(&mut self, flags: u16, mapq: u8, ref_id: i32, mate_ref_id: i32) {
        let qc_fail = (flags & FLAG_QCFAIL) != 0;
        let c = if qc_fail {
            &mut self.fail
        } else {
            &mut self.pass
        };
        update_counts(c, flags, mapq, ref_id, mate_ref_id);
    }
}

fn update_counts(c: &mut Counts, flags: u16, mapq: u8, ref_id: i32, mate_ref_id: i32) {
    c.total += 1;
    let is_secondary = (flags & FLAG_SECONDARY) != 0;
    let is_supp = (flags & FLAG_SUPPLEMENTARY) != 0;
    let is_primary = !is_secondary && !is_supp;
    if is_secondary {
        c.secondary += 1;
    }
    if is_supp {
        c.supplementary += 1;
    }
    if is_primary {
        c.primary += 1;
    }
    let dup = (flags & FLAG_DUP) != 0;
    if dup {
        c.duplicates += 1;
        if is_primary {
            c.primary_duplicates += 1;
        }
    }
    let unmapped = (flags & FLAG_UNMAP) != 0;
    if !unmapped {
        c.mapped += 1;
        if is_primary {
            c.primary_mapped += 1;
        }
    }
    if (flags & FLAG_PAIRED) != 0 {
        c.paired += 1;
        if (flags & FLAG_READ1) != 0 {
            c.read1 += 1;
        }
        if (flags & FLAG_READ2) != 0 {
            c.read2 += 1;
        }
        if !unmapped {
            if (flags & FLAG_MUNMAP) == 0 {
                c.with_itself_and_mate_mapped += 1;
                if (flags & FLAG_PROPER_PAIR) != 0 {
                    c.properly_paired += 1;
                }
                // samtools applies `mtid != tid` literally even when mtid is -1.
                if ref_id != mate_ref_id {
                    c.mate_diff_chr += 1;
                    if mapq >= 5 {
                        c.mate_diff_chr_mq5 += 1;
                    }
                }
            } else {
                c.singletons += 1;
            }
        }
    }
}

pub fn run(args: FlagstatArgs) -> Result<()> {
    install_thread_pool(args.threads);
    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let mut totals = FlagstatTotals::default();
    let mut rec = RecordSummary::default();
    while reader.read_summary(&mut rec)? {
        totals.update(
            rec.flags,
            rec.mapping_quality.unwrap_or(0),
            rec.tid(),
            rec.mate_tid(),
        );
    }
    match args.output_fmt.as_str() {
        "json" => print_flagstat_json(&totals),
        _ => print_flagstat(&totals),
    }
    Ok(())
}

/// samtools' `-O json` schema: two objects keyed `"QC-passed reads"` and
/// `"QC-failed reads"`, each holding the same field names the plain-text report
/// uses, with the percentages as their own `… %` entries. The old flat
/// `{key: {pass, fail}}` shape parsed as valid JSON but no samtools consumer
/// could read it.
fn print_flagstat_json(t: &FlagstatTotals) {
    println!("{{");
    println!("{},", counts_json("QC-passed reads", &t.pass));
    println!("{}", counts_json("QC-failed reads", &t.fail));
    println!("}}");
}

fn counts_json(label: &str, c: &Counts) -> String {
    /// Percentages are `null` where samtools prints `N/A`.
    fn ratio(num: u64, den: u64) -> String {
        if den == 0 {
            "null".to_string()
        } else {
            format!("{:.6}", (num as f64) * 100.0 / (den as f64))
        }
    }
    let fields = [
        format!("\"total\": {}", c.total),
        format!("\"primary\": {}", c.primary),
        format!("\"secondary\": {}", c.secondary),
        format!("\"supplementary\": {}", c.supplementary),
        format!("\"duplicates\": {}", c.duplicates),
        format!("\"primary duplicates\": {}", c.primary_duplicates),
        format!("\"mapped\": {}", c.mapped),
        format!("\"mapped %\": {}", ratio(c.mapped, c.total)),
        format!("\"primary mapped\": {}", c.primary_mapped),
        format!(
            "\"primary mapped %\": {}",
            ratio(c.primary_mapped, c.primary)
        ),
        format!("\"paired in sequencing\": {}", c.paired),
        format!("\"read1\": {}", c.read1),
        format!("\"read2\": {}", c.read2),
        format!("\"properly paired\": {}", c.properly_paired),
        format!(
            "\"properly paired %\": {}",
            ratio(c.properly_paired, c.paired)
        ),
        format!(
            "\"with itself and mate mapped\": {}",
            c.with_itself_and_mate_mapped
        ),
        format!("\"singletons\": {}", c.singletons),
        format!("\"singletons %\": {}", ratio(c.singletons, c.paired)),
        format!(
            "\"with mate mapped to a different chr\": {}",
            c.mate_diff_chr
        ),
        format!(
            "\"with mate mapped to a different chr (mapQ>=5)\": {}",
            c.mate_diff_chr_mq5
        ),
    ];
    let body = fields
        .iter()
        .map(|f| format!("    {f}"))
        .collect::<Vec<_>>()
        .join(",\n");
    format!("  \"{label}\": {{\n{body}\n  }}")
}

fn pct(num: u64, den: u64) -> String {
    if den == 0 {
        "N/A".to_string()
    } else {
        format!("{:.2}%", (num as f64) * 100.0 / (den as f64))
    }
}

fn print_flagstat(t: &FlagstatTotals) {
    let p = &t.pass;
    let f = &t.fail;
    println!(
        "{} + {} in total (QC-passed reads + QC-failed reads)",
        p.total, f.total
    );
    println!("{} + {} primary", p.primary, f.primary);
    println!("{} + {} secondary", p.secondary, f.secondary);
    println!("{} + {} supplementary", p.supplementary, f.supplementary);
    println!("{} + {} duplicates", p.duplicates, f.duplicates);
    println!(
        "{} + {} primary duplicates",
        p.primary_duplicates, f.primary_duplicates
    );
    println!(
        "{} + {} mapped ({} : {})",
        p.mapped,
        f.mapped,
        pct(p.mapped, p.total),
        pct(f.mapped, f.total)
    );
    println!(
        "{} + {} primary mapped ({} : {})",
        p.primary_mapped,
        f.primary_mapped,
        pct(p.primary_mapped, p.primary),
        pct(f.primary_mapped, f.primary)
    );
    println!("{} + {} paired in sequencing", p.paired, f.paired);
    println!("{} + {} read1", p.read1, f.read1);
    println!("{} + {} read2", p.read2, f.read2);
    println!(
        "{} + {} properly paired ({} : {})",
        p.properly_paired,
        f.properly_paired,
        pct(p.properly_paired, p.paired),
        pct(f.properly_paired, f.paired)
    );
    println!(
        "{} + {} with itself and mate mapped",
        p.with_itself_and_mate_mapped, f.with_itself_and_mate_mapped
    );
    println!(
        "{} + {} singletons ({} : {})",
        p.singletons,
        f.singletons,
        pct(p.singletons, p.paired),
        pct(f.singletons, f.paired)
    );
    println!(
        "{} + {} with mate mapped to a different chr",
        p.mate_diff_chr, f.mate_diff_chr
    );
    println!(
        "{} + {} with mate mapped to a different chr (mapQ>=5)",
        p.mate_diff_chr_mq5, f.mate_diff_chr_mq5
    );
}
