use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as AlignmentWrite;

use crate::cli::HeadArgs;
use crate::io::BamReader;

pub fn run(args: HeadArgs) -> Result<()> {
    let mut reader = BamReader::open(&args.input).context("open input")?;
    let header = reader.header().clone();

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    // Emit header.
    {
        let mut sam_writer = noodles_sam::io::Writer::new(&mut out);
        sam_writer.write_header(&header)?;
    }

    if args.n == 0 {
        return Ok(());
    }
    let mut sam_writer = noodles_sam::io::Writer::new(&mut out);
    let mut rec = RecordBuf::default();
    let mut n: u64 = 0;
    while reader.read_record_buf(&mut rec)? {
        sam_writer.write_alignment_record(&header, &rec)?;
        n += 1;
        if n >= args.n {
            break;
        }
    }
    Ok(())
}
