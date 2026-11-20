use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, ensure};

pub(crate) struct FastqRecord<'a> {
    pub name: &'a str,
    pub comment: &'a str,
    pub sequence: &'a [u8],
    pub quality: &'a [u8],
    pub sequence_offset: u64,
    pub quality_offset: u64,
}

pub(crate) struct FastqInput {
    input: BufReader<File>,
    path: PathBuf,
    record_number: u64,
    offset: u64,
    lines: [String; 4],
}

impl FastqInput {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open FASTQ {}", path.display()))?;
        Ok(Self {
            input: BufReader::with_capacity(1 << 20, file),
            path: path.to_owned(),
            record_number: 0,
            offset: 0,
            lines: std::array::from_fn(|_| String::new()),
        })
    }

    pub(crate) fn next(&mut self) -> Result<Option<FastqRecord<'_>>> {
        let start = self.offset;
        self.record_number += 1;
        if !self.read_lines().with_context(|| {
            format!(
                "{}: FASTQ record {}",
                self.path.display(),
                self.record_number
            )
        })? {
            return Ok(None);
        }
        self.parse(start).map(Some).with_context(|| {
            format!(
                "{}: FASTQ record {}",
                self.path.display(),
                self.record_number
            )
        })
    }

    fn read_lines(&mut self) -> Result<bool> {
        for (index, line) in self.lines.iter_mut().enumerate() {
            line.clear();
            let len = self.input.read_line(line)?;
            if len == 0 {
                ensure!(index == 0, "truncated FASTQ record");
                return Ok(false);
            }
            self.offset += len as u64;
        }
        Ok(true)
    }

    fn parse(&self, start: u64) -> Result<FastqRecord<'_>> {
        let [header, sequence, plus, quality] = &self.lines;
        let sequence_offset = start + header.len() as u64;
        let quality_offset = sequence_offset + sequence.len() as u64 + plus.len() as u64;
        let header = header
            .trim_end_matches(['\r', '\n'])
            .strip_prefix('@')
            .context("missing @ header marker")?;
        let (name, comment) = header
            .split_once(char::is_whitespace)
            .unwrap_or((header, ""));
        ensure!(!name.is_empty(), "missing read name");
        ensure!(plus.starts_with('+'), "missing + separator");
        let sequence = sequence.trim_end_matches(['\r', '\n']).as_bytes();
        let quality = quality.trim_end_matches(['\r', '\n']).as_bytes();
        ensure!(
            sequence.len() == quality.len(),
            "sequence and quality lengths differ"
        );
        ensure!(
            quality.iter().all(|q| (33..=126).contains(q)),
            "invalid FASTQ quality"
        );
        Ok(FastqRecord {
            name,
            comment,
            sequence,
            quality,
            sequence_offset,
            quality_offset,
        })
    }
}
