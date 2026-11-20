use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use tempfile::{NamedTempFile, TempPath};

pub(super) struct PendingOutput {
    staging: TempPath,
    destination: PathBuf,
}

impl PendingOutput {
    pub(super) fn create(path: &Path) -> Result<(File, Option<Self>)> {
        let metadata = match fs::metadata(path) {
            Ok(metadata) => Some(metadata),
            Err(err) if err.kind() == io::ErrorKind::NotFound => None,
            Err(err) => return Err(err).with_context(|| format!("inspect {}", path.display())),
        };
        if metadata.as_ref().is_some_and(|m| !m.is_file()) {
            let file = File::options()
                .write(true)
                .open(path)
                .with_context(|| format!("open output stream {}", path.display()))?;
            return Ok((file, None));
        }
        let destination = if fs::symlink_metadata(path).is_ok_and(|m| m.is_symlink()) {
            fs::canonicalize(path)
                .with_context(|| format!("resolve output link {}", path.display()))?
        } else {
            path.to_owned()
        };
        let (file, staging) = NamedTempFile::new_in(output_parent(&destination))
            .with_context(|| format!("create temporary output for {}", path.display()))?
            .into_parts();
        Ok((
            file,
            Some(Self {
                staging,
                destination,
            }),
        ))
    }

    pub(super) fn publish(self) -> Result<()> {
        self.staging
            .persist(&self.destination)
            .with_context(|| format!("publish {}", self.destination.display()))
    }
}

pub(crate) struct BufferedOutput {
    writer: BufWriter<Box<dyn Write + Send>>,
    pending: Option<PendingOutput>,
}

impl BufferedOutput {
    pub(crate) fn create(path: Option<&Path>) -> Result<Self> {
        let (sink, pending): (Box<dyn Write + Send>, _) =
            if let Some(path) = path.filter(|p| *p != Path::new("-")) {
                let (file, pending) = PendingOutput::create(path)?;
                (Box::new(file), pending)
            } else {
                (Box::new(io::stdout()), None)
            };
        Ok(Self {
            writer: BufWriter::with_capacity(1 << 16, sink),
            pending,
        })
    }

    pub(crate) fn finish(mut self) -> Result<()> {
        self.writer.flush().context("flush output")?;
        drop(self.writer);
        if let Some(pending) = self.pending {
            pending.publish()?;
        }
        Ok(())
    }
}

impl Write for BufferedOutput {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.writer.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

fn output_parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

/// Compare destinations even when the file has yet to be created.
pub(crate) fn output_key(path: &Path) -> Result<PathBuf> {
    let resolved = match fs::canonicalize(path) {
        Ok(path) => path,
        Err(err) if err.kind() == io::ErrorKind::NotFound => {
            let name = path.file_name().context("output path needs a file name")?;
            fs::canonicalize(output_parent(path))?.join(name)
        }
        Err(err) => return Err(err).with_context(|| format!("resolve {}", path.display())),
    };
    #[cfg(windows)]
    let resolved = PathBuf::from(resolved.to_string_lossy().to_lowercase());
    Ok(resolved)
}

pub(crate) fn ensure_distinct_output(input: &Path, output: &Path) -> Result<PathBuf> {
    let key = output_key(output)?;
    let same_file = match same_file::is_same_file(input, output) {
        Ok(same) => same,
        Err(err) if err.kind() == io::ErrorKind::NotFound => false,
        Err(err) => return Err(err).context("compare input and output files"),
    };
    anyhow::ensure!(
        !same_file && output_key(input)? != key,
        "output {} refers to input {}",
        output.display(),
        input.display()
    );
    Ok(key)
}
