use std::fs::{self, File};
use std::io;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use tempfile::{NamedTempFile, TempPath};

pub(super) struct PendingOutput {
    staging: TempPath,
    destination: PathBuf,
}

impl PendingOutput {
    pub(super) fn create(path: &Path) -> Result<(File, Self)> {
        let directory = output_parent(path);
        let (file, staging) = NamedTempFile::new_in(directory)
            .with_context(|| format!("create temporary output for {}", path.display()))?
            .into_parts();
        Ok((
            file,
            Self {
                staging,
                destination: path.to_owned(),
            },
        ))
    }

    pub(super) fn publish(self) -> Result<()> {
        self.staging
            .persist(&self.destination)
            .with_context(|| format!("publish {}", self.destination.display()))
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
