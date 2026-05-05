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

pub fn record_start_0(record: &bam::Record) -> Option<u64> {
    record.alignment_start()?.ok().map(|p| p.get() as u64 - 1)
}

pub fn record_is_owned_by_window(record: &bam::Record, window: &BamWindow) -> bool {
    let Some(start) = record_start_0(record) else {
        return false;
    };
    start >= window.start as u64 && start < window.end as u64
}

pub fn scan_bam_windows_full<T, Make, Visit, Finalize, Merge>(
    bam_path: &str,
    bai_path: &str,
    chunk_size: usize,
    show_progress: bool,
    make_state: Make,
    visit: Visit,
    finalize: Finalize,
    merge: Merge,
) -> Result<T>
where
    T: Send,
    Make: Fn() -> T + Sync + Send,
    Visit: Fn(&mut T, &sam::Header, &BamWindow, &bam::Record) + Sync + Send,
    Finalize: Fn(&mut T, &BamWindow) + Sync + Send,
    Merge: Fn(T, T) -> T + Sync + Send,
{
    let header = {
        let file = File::open(bam_path).with_context(|| format!("could not open {bam_path}"))?;
        let mut reader = bam::io::Reader::new(file);
        reader.read_header()?
    };
    let windows = generate_windows(&header, chunk_size);
    let bai = bam::bai::read(bai_path)?;

    let pb = if show_progress {
        let pb = ProgressBar::new(windows.len() as u64);
        pb.set_style(ProgressStyle::default_bar().template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} windows ({eta})",
        )?);
        Some(pb)
    } else {
        None
    };

    let states = windows
        .par_iter()
        .map(|window| -> Result<T> {
            let file =
                File::open(bam_path).with_context(|| format!("could not open {bam_path}"))?;
            let mut reader = bam::io::indexed_reader::Builder::default()
                .set_index(bai.clone())
                .build_from_reader(file)?;
            let _ = reader.read_header()?;
            let region: noodles::core::Region =
                format!("{}:{}-{}", window.chrom, window.start + 1, window.end).parse()?;
            let mut state = make_state();
            let query = reader.query(&header, &region)?;
            for result in query {
                let record = result?;
                if record_is_owned_by_window(&record, window) {
                    visit(&mut state, &header, window, &record);
                }
            }
            finalize(&mut state, window);
            if let Some(pb) = &pb {
                pb.inc(1);
            }
            Ok(state)
        })
        .collect::<Result<Vec<_>>>()?;

    if let Some(pb) = pb {
        pb.finish_with_message("Window scan finished.");
    }

    Ok(states.into_par_iter().reduce(make_state, merge))
}

pub fn scan_bam_windows<T, Make, Visit, Merge>(
    bam_path: &str,
    bai_path: &str,
    chunk_size: usize,
    show_progress: bool,
    make_state: Make,
    visit: Visit,
    merge: Merge,
) -> Result<T>
where
    T: Send,
    Make: Fn() -> T + Sync + Send,
    Visit: Fn(&mut T, &sam::Header, &BamWindow, &bam::Record) + Sync + Send,
    Merge: Fn(T, T) -> T + Sync + Send,
{
    scan_bam_windows_full(
        bam_path,
        bai_path,
        chunk_size,
        show_progress,
        make_state,
        visit,
        |_state, _window| {},
        merge,
    )
}
