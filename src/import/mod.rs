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
    let mut writer = BamWriter::create(Some(&args.output), header, fmt)?;
    writer.write_header()?;

    let mut r1 = FastqIter::open(&args.r1)?;
    let mut r2 = match &args.r2 {
        Some(p) => Some(FastqIter::open(p)?),
        None => None,
    };

    let rg_id = args.rg_id.clone();
    while let Some(rec1) = r1.next()? {
        let flag_r1: u16 = if r2.is_some() {
            0x1 | 0x40 | 0x4 | 0x8
        } else {
            0x4
        };
        emit(&mut writer, &rec1, flag_r1, &rg_id, args.casava, 1)?;
        if let Some(r2) = r2.as_mut()
            && let Some(rec2) = r2.next()?
        {
            let flag_r2 = 0x1 | 0x80 | 0x4 | 0x8;
            emit(&mut writer, &rec2, flag_r2, &rg_id, args.casava, 2)?;
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
}

impl FastqIter {
    fn open(p: &std::path::Path) -> Result<Self> {
        Ok(Self {
            reader: BufReader::with_capacity(1 << 20, File::open(p)?),
        })
    }
    fn next(&mut self) -> Result<Option<FastqRec>> {
        let mut header = String::new();
        if self.reader.read_line(&mut header)? == 0 {
            return Ok(None);
        }
        let mut seq = String::new();
        if self.reader.read_line(&mut seq)? == 0 {
            return Ok(None);
        }
        let mut plus = String::new();
        if self.reader.read_line(&mut plus)? == 0 {
            return Ok(None);
        }
        let mut qual = String::new();
        if self.reader.read_line(&mut qual)? == 0 {
            return Ok(None);
        }
        let header = header.trim_end_matches(['\n', '\r']);
        let seq = seq.trim_end_matches(['\n', '\r']);
        let qual = qual.trim_end_matches(['\n', '\r']);
        let stripped = header.strip_prefix('@').unwrap_or(header);
        let mut parts = stripped.splitn(2, char::is_whitespace);
        let name = parts.next().unwrap_or("").to_string();
        // Strip /1 /2 suffix.
        let name = name
            .trim_end_matches("/1")
            .trim_end_matches("/2")
            .to_string();
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
            name: name.into_bytes(),
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
    _mate_num: u8,
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
