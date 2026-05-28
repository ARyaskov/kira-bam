use std::path::Path;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OutputFormat {
    Sam,
    Bam,
    UncompressedBam,
    Cram,
}

impl OutputFormat {
    pub fn from_path_ext(path: &Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some("sam") => OutputFormat::Sam,
            Some("bam") => OutputFormat::Bam,
            Some("cram") => OutputFormat::Cram,
            _ => OutputFormat::Bam,
        }
    }
}
