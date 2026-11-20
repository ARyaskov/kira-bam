use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as AlignmentWrite;

use crate::cli::HeadArgs;
use crate::io::BamReader;
use crate::io::BufferedOutput;

pub fn run(args: HeadArgs) -> Result<()> {
    let mut reader = BamReader::open_with_reference(&args.input, args.reference.as_deref())
        .context("open input")?;
    let header = reader.header().clone();

    let mut out = BufferedOutput::create(args.output.as_deref())?;

    // Emit header.
    {
        let mut sam_writer = noodles_sam::io::Writer::new(&mut out);
        sam_writer.write_header(&header)?;
    }

    if args.n == 0 {
        return out.finish();
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
    out.finish()
}
