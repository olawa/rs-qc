use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use noodles::bam;
use noodles::sam;
use rayon::prelude::*;
use std::fs::File;
use std::num::NonZeroUsize;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct BamScanConfig {
    pub threads: usize,
    pub show_progress: bool,
}

#[derive(Clone, Debug)]
pub struct BamWindow {
    pub chrom: String,
    pub chrom_norm: String,
    pub start: u32,
    pub end: u32,
}

pub fn generate_windows(header: &sam::Header, chunk_size: usize) -> Vec<BamWindow> {
    let mut windows = Vec::new();
    for (name, seq) in header.reference_sequences() {
        let chrom = String::from_utf8_lossy(name.as_ref()).to_string();
        let chrom_norm = crate::models::normalize_chrom(&chrom).into_owned();
        let len = seq.length().get() as u32;
        let mut start = 0;
        while start < len {
            let end = (start + chunk_size as u32).min(len);
            windows.push(BamWindow {
                chrom: chrom.clone(),
                chrom_norm: chrom_norm.clone(),
                start,
                end,
            });
            start = end;
        }
    }
    windows
}

pub fn generate_windows_for_range(
    chrom: &str,
    chrom_norm: &str,
    start: u64,
    end: u64,
    window_size: usize,
) -> Vec<BamWindow> {
    let mut windows = Vec::new();
    let mut curr = start as u32;
    let end_u32 = end as u32;
    while curr < end_u32 {
        let win_end = (curr + window_size as u32).min(end_u32);
        windows.push(BamWindow {
            chrom: chrom.to_string(),
            chrom_norm: chrom_norm.to_string(),
            start: curr,
            end: win_end,
        });
        curr = win_end;
    }
    windows
}

pub fn find_bai_path(bam_path: &str) -> Option<String> {
    if Path::new(&format!("{bam_path}.bai")).exists() {
        Some(format!("{bam_path}.bai"))
    } else {
        bam_path
            .strip_suffix(".bam")
            .map(|s| format!("{s}.bai"))
            .filter(|s| Path::new(s).exists())
    }
}

pub fn scan_bam_stream<T, F>(
    bam_path: &str,
    config: &BamScanConfig,
    mut state: T,
    mut visit: F,
) -> Result<T>
where
    F: FnMut(&mut T, &sam::Header, &bam::Record),
{
    let file = File::open(bam_path).with_context(|| format!("could not open {bam_path}"))?;
    let worker_count = NonZeroUsize::new(config.threads).unwrap_or(NonZeroUsize::MIN);
    let bgzf_reader = noodles::bgzf::MultithreadedReader::with_worker_count(worker_count, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);
    let header = reader.read_header()?;

    let pb = if config.show_progress {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {pos} records scanned")?,
        );
        Some(pb)
    } else {
        None
    };

    for result in reader.records() {
        let record = result?;
        visit(&mut state, &header, &record);
        if let Some(pb) = &pb {
            pb.inc(1);
        }
    }

    if let Some(pb) = pb {
        pb.finish_with_message("BAM scan finished.");
    }

    Ok(state)
}

