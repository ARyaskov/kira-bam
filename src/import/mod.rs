use anyhow::{Context, Result};
use bstr::BString;
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::header::record::value::{Map, map::ReadGroup, map::read_group::tag as rg_tag};

use crate::cli::ImportArgs;
use crate::fastq_input::{FastqInput, FastqRecord};
use crate::io::{BamWriter, PgInfo, append_pg};
use crate::types::OutputFormat;

pub fn run(args: ImportArgs) -> Result<()> {
    let mut header = Header::default();
    let mut rg = Map::<ReadGroup>::builder();
    if let Some(s) = &args.sample {
        rg = rg.insert(rg_tag::SAMPLE, s.as_bytes());
    }
    let rg = rg.build().context("build @RG")?;
    header
        .read_groups_mut()
        .insert(BString::from(args.rg_id.as_bytes()), rg);
    append_pg(&mut header, &PgInfo::new("import", !args.no_pg))?;

    let fmt = if args.uncompressed {
        OutputFormat::UncompressedBam
    } else {
        OutputFormat::Bam
    };
    let mut r1 = FastqInput::open(&args.r1)?;
    let mut r2 = match &args.r2 {
        Some(p) => Some(FastqInput::open(p)?),
        None => None,
    };
    let mut writer = BamWriter::create(Some(&args.output), header, fmt)?;
    writer.write_header()?;

    let rg_id = args.rg_id.clone();
    if let Some(r2) = &mut r2 {
        for pair in 1u64.. {
            match (r1.next()?, r2.next()?) {
                (None, None) => break,
                (Some(rec1), Some(rec2)) => {
                    anyhow::ensure!(
                        template_name(rec1.name) == template_name(rec2.name),
                        "FASTQ pair {pair}: mate names differ ({:?}, {:?})",
                        rec1.name,
                        rec2.name
                    );
                    emit(&mut writer, &rec1, 0x4d, &rg_id, args.casava)?;
                    emit(&mut writer, &rec2, 0x8d, &rg_id, args.casava)?;
                }
                _ => anyhow::bail!("FASTQ pair {pair}: mate files have different record counts"),
            }
        }
    } else {
        while let Some(rec) = r1.next()? {
            emit(&mut writer, &rec, 0x4, &rg_id, args.casava)?;
        }
    }
    writer.finish()?;
    Ok(())
}

fn template_name(name: &str) -> &str {
    name.strip_suffix("/1")
        .or_else(|| name.strip_suffix("/2"))
        .unwrap_or(name)
}

fn emit(
    writer: &mut BamWriter,
    rec: &FastqRecord<'_>,
    flag_bits: u16,
    rg_id: &str,
    casava: bool,
) -> Result<()> {
    let mut rb = RecordBuf::default();
    let name = template_name(rec.name);
    anyhow::ensure!(!name.is_empty(), "missing read name");
    *rb.name_mut() = Some(BString::from(name.as_bytes()));
    *rb.flags_mut() = Flags::from(flag_bits);
    *rb.sequence_mut() = rec.sequence.to_vec().into();
    let q: Vec<u8> = rec.quality.iter().map(|&c| c - 33).collect();
    *rb.quality_scores_mut() = q.into();
    rb.data_mut().insert(
        Tag::READ_GROUP,
        Value::String(BString::from(rg_id.as_bytes())),
    );
    let barcode = rec.comment.rsplit(':').next().unwrap_or("");
    if casava
        && !barcode.is_empty()
        && barcode
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'+' | b'-'))
    {
        rb.data_mut().insert(
            Tag::new(b'B', b'C'),
            Value::String(BString::from(barcode.as_bytes())),
        );
    }
    writer.write_record_buf(&rb)?;
    Ok(())
}
