use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Result};
use bstr::BString;
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::header::record::value::{Map, map::ReadGroup, map::read_group::tag as rg_tag};

use crate::cli::ImportArgs;
use crate::io::{BamWriter, PgInfo, append_pg, install_thread_pool};
use crate::types::OutputFormat;

pub fn run(args: ImportArgs) -> Result<()> {
    install_thread_pool(args.threads);

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
    let mut r1 = FastqIter::open(&args.r1)?;
    let mut r2 = match &args.r2 {
        Some(p) => Some(FastqIter::open(p)?),
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
                        rec1.name == rec2.name,
                        "FASTQ pair {pair}: mate names differ ({:?}, {:?})",
                        String::from_utf8_lossy(&rec1.name),
                        String::from_utf8_lossy(&rec2.name)
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

struct FastqRec {
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    barcode: Option<Vec<u8>>,
}

struct FastqIter {
    reader: BufReader<File>,
    path: std::path::PathBuf,
    record: u64,
}

impl FastqIter {
    fn open(p: &std::path::Path) -> Result<Self> {
        Ok(Self {
            reader: BufReader::with_capacity(
                1 << 20,
                File::open(p).with_context(|| format!("open FASTQ {}", p.display()))?,
            ),
            path: p.to_owned(),
            record: 0,
        })
    }
    fn next(&mut self) -> Result<Option<FastqRec>> {
        self.record += 1;
        let record = self.record;
        self.read_record()
            .with_context(|| format!("{}: FASTQ record {record}", self.path.display()))
    }

    fn read_record(&mut self) -> Result<Option<FastqRec>> {
        let mut header = String::new();
        if self.reader.read_line(&mut header)? == 0 {
            return Ok(None);
        }
        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();
        anyhow::ensure!(
            self.reader.read_line(&mut seq)? != 0
                && self.reader.read_line(&mut plus)? != 0
                && self.reader.read_line(&mut qual)? != 0,
            "truncated FASTQ record"
        );
        let header = header.trim_end_matches(['\n', '\r']);
        let seq = seq.trim_end_matches(['\n', '\r']);
        let qual = qual.trim_end_matches(['\n', '\r']);
        let stripped = header
            .strip_prefix('@')
            .context("missing @ header marker")?;
        anyhow::ensure!(plus.starts_with('+'), "missing + separator");
        anyhow::ensure!(
            seq.len() == qual.len(),
            "sequence and quality lengths differ"
        );
        anyhow::ensure!(
            qual.bytes().all(|q| (33..=126).contains(&q)),
            "invalid FASTQ quality"
        );
        let mut parts = stripped.splitn(2, char::is_whitespace);
        let name = parts.next().unwrap_or("");
        let name = name
            .strip_suffix("/1")
            .or_else(|| name.strip_suffix("/2"))
            .unwrap_or(name);
        anyhow::ensure!(!name.is_empty(), "missing read name");
        let rest = parts.next().unwrap_or("");
        // CASAVA: 1:N:0:BARCODE → pull BC.
        let mut barcode: Option<Vec<u8>> = None;
        if !rest.is_empty() {
            let segs: Vec<&str> = rest.split(':').collect();
            if let Some(bc) = segs.last()
                && !bc.is_empty()
                && bc
                    .chars()
                    .all(|c| c.is_ascii_alphanumeric() || c == '+' || c == '-')
            {
                barcode = Some(bc.as_bytes().to_vec());
            }
        }
        Ok(Some(FastqRec {
            name: name.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: qual.as_bytes().to_vec(),
            barcode,
        }))
    }
}

fn emit(
    writer: &mut BamWriter,
    rec: &FastqRec,
    flag_bits: u16,
    rg_id: &str,
    casava: bool,
) -> Result<()> {
    let mut rb = RecordBuf::default();
    *rb.name_mut() = Some(BString::from(rec.name.clone()));
    *rb.flags_mut() = Flags::from(flag_bits);
    *rb.sequence_mut() = rec.seq.clone().into();
    let q: Vec<u8> = rec.qual.iter().map(|&c| c.saturating_sub(33)).collect();
    *rb.quality_scores_mut() = q.into();
    rb.data_mut().insert(
        Tag::READ_GROUP,
        Value::String(BString::from(rg_id.as_bytes())),
    );
    if casava && let Some(bc) = &rec.barcode {
        rb.data_mut().insert(
            Tag::new(b'B', b'C'),
            Value::String(BString::from(bc.clone())),
        );
    }
    writer.write_record_buf(&rb)?;
    Ok(())
}
