use std::path::PathBuf;

use clap::Parser;

/// Text alignment viewer. Non-interactive: prints a fixed-width pile to stdout
/// (samtools `tview -d T` style). For an actual TUI we'd need a terminal lib —
/// out of scope while we target Windows + Unix without an extra dep.
#[derive(Parser, Debug)]
pub struct TviewArgs {
    /// Indexed BAM file.
    #[arg(value_name = "IN")]
    pub input: PathBuf,

    /// Reference FASTA (with .fai); optional.
    #[arg(value_name = "FASTA")]
    pub reference: Option<PathBuf>,

    /// Position to display, `chr:pos` (1-based).
    #[arg(short = 'p', long = "position", required = true)]
    pub position: String,

    /// Display width in bases (default 80).
    #[arg(short = 'w', long = "width", default_value_t = 80)]
    pub width: usize,

    /// Max number of reads to stack vertically.
    #[arg(short = 'd', long = "depth", default_value_t = 30)]
    pub depth: usize,

    /// Suppress ANSI colour escapes.
    #[arg(long = "no-color")]
    pub no_color: bool,
}
