use super::config::RnaQcConfig;
use super::state::RnaWorkerState;
use crate::analysis::bam_scan::{find_bai_path, generate_windows};
use crate::analysis::contamination::ContaminantIndex;
use crate::analysis::feature_index::{ChromFeatureIndex, FeatureIndex};
use crate::analysis::index::{AnnotationIndex, DenseMap, Hits};
use crate::analysis::types::AnalysisType;
use crate::io::bam::{alignment_start_0, for_each_aligned_block, match_span, reference_span};
use crate::models::{is_autosomal_chrom, normalize_chrom};
use anyhow::{Context, Result};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use noodles::bam;
use noodles::sam;
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::HashSet;
use std::fs::File;
use std::sync::Arc;

thread_local! {
    static BAM_READER: RefCell<Option<(String, bam::io::IndexedReader<noodles::bgzf::Reader<File>>)>> = RefCell::new(None);
}

pub(crate) fn normalized_name_set(names: &[String]) -> HashSet<String> {
    names
        .iter()
        .map(|n| normalize_chrom(n).into_owned())
        .collect()
}

pub(crate) fn autosomal_coverage_index(index: &AnnotationIndex) -> AnnotationIndex {
    let genes = index
        .genes
        .iter()
        .filter(|gene| is_autosomal_chrom(&gene.chrom))
        .cloned()
        .collect();
    AnnotationIndex::new(genes, false)
}

pub(crate) fn normalize_thread_count(threads: usize) -> usize {
    threads.max(1)
}

pub(crate) fn window_owns_record_start(pos: u64, win_start: u64, win_end: u64) -> bool {
    pos >= win_start && pos < win_end
}

pub(crate) fn is_mtdna_chrom_norm(chrom: &str) -> bool {
    chrom == "m" || chrom == "mt" || chrom == "chrm" || chrom == "chrmt"
}

pub(crate) fn is_rdna_chrom_norm(chrom: &str, rdna_contigs: &HashSet<String>) -> bool {
    rdna_contigs.contains(chrom)
}

pub(crate) struct ScanContext<'a> {
    pub(crate) bam_path: &'a str,
    pub(crate) header: &'a sam::Header,
    pub(crate) dense_maps: &'a std::collections::HashMap<String, Arc<DenseMap>>,
    pub(crate) coverage_index: &'a AnnotationIndex,
    pub(crate) feature_index: Option<&'a FeatureIndex>,
    pub(crate) rdna_contigs: &'a HashSet<String>,
    pub(crate) rdna_intervals: Option<&'a ContaminantIndex>,
    pub(crate) config: &'a RnaQcConfig,
    pub(crate) threads: usize,
}

pub(crate) fn scan_bam_coverage(ctx: &ScanContext<'_>) -> Result<RnaWorkerState> {
    if let Some(bai_path) = find_bai_path(ctx.bam_path) {
        println!("  - BAI Index found. Using high-performance parallel dense scan.");
        let windows = generate_windows(ctx.header, 10_000_000);
        let pb = ProgressBar::new(windows.len() as u64);
        pb.set_style(ProgressStyle::default_bar().template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Windows ({eta})",
        )?);
        let res = scan_windows_with_pb(ctx, &bai_path, &windows, Some(pb.clone()));
        pb.finish_and_clear();
        res
    } else {
        println!("  - No BAI Index found. Falling back to sequential single-threaded scan.");
        scan_sequential_bam(ctx)
    }
}

pub(crate) fn scan_windows(
    ctx: &ScanContext<'_>,
    bai_path: &str,
    windows: &[crate::analysis::bam_scan::BamWindow],
) -> Result<RnaWorkerState> {
    scan_indexed_bam(ctx, bai_path, windows, None)
}

pub(crate) fn scan_windows_with_pb(
    ctx: &ScanContext<'_>,
    bai_path: &str,
    windows: &[crate::analysis::bam_scan::BamWindow],
    pb: Option<ProgressBar>,
) -> Result<RnaWorkerState> {
    scan_indexed_bam(ctx, bai_path, windows, pb)
}

fn scan_indexed_bam(
    ctx: &ScanContext<'_>,
    bai_path: &str,
    windows: &[crate::analysis::bam_scan::BamWindow],
    pb: Option<ProgressBar>,
) -> Result<RnaWorkerState> {
    let bai = bam::bai::read(bai_path)?;

    let bai_arc = Arc::new(bai);
    let qc_sample_size = ctx.config.qc_sample_size.div_ceil(ctx.threads);
    let empty_state = || RnaWorkerState::new(qc_sample_size);

    windows
        .par_iter()
        .fold(
            || Ok(empty_state()),
            |state_res, window| -> Result<RnaWorkerState> {
                let mut state = state_res?;
                let window_res = BAM_READER.with(|cell| -> Result<()> {
                    let mut opt = cell.borrow_mut();
                    let needs_reader = opt
                        .as_ref()
                        .map(|(path, _)| path != ctx.bam_path)
                        .unwrap_or(true);
                    if needs_reader {
                        let f = File::open(ctx.bam_path)
                            .with_context(|| format!("failed to open BAM {}", ctx.bam_path))?;
                        let mut reader = bam::io::indexed_reader::Builder::default()
                            .set_index(bai_arc.as_ref().clone())
                            .build_from_reader(f)
                            .with_context(|| {
                                format!("failed to build indexed reader for {}", ctx.bam_path)
                            })?;
                        reader.read_header().with_context(|| {
                            format!("failed to read BAM header for {}", ctx.bam_path)
                        })?;
                        *opt = Some((ctx.bam_path.to_string(), reader));
                    }

                    let reader = &mut opt.as_mut().unwrap().1;
                    scan_indexed_window(ctx, reader, window, &mut state)
                });
                window_res?;
                if let Some(pb) = &pb {
                    pb.inc(1);
                }
                Ok(state)
            },
        )
        .reduce(
            || Ok(empty_state()),
            |a, b| match (a, b) {
                (Ok(left), Ok(right)) => Ok(left.merge(right)),
                (Err(e), _) => Err(e),
                (_, Err(e)) => Err(e),
            },
        )
}

fn scan_indexed_window(
    ctx: &ScanContext<'_>,
    reader: &mut bam::io::IndexedReader<noodles::bgzf::Reader<File>>,
    window: &crate::analysis::bam_scan::BamWindow,
    state: &mut RnaWorkerState,
) -> Result<()> {
    let chrom = &window.chrom;
    let chrom_norm = &window.chrom_norm;
    let maybe_dense = ctx.dense_maps.get(chrom_norm);
    let maybe_feature = ctx
        .feature_index
        .and_then(|index| index.chroms.get(chrom_norm));
    let is_mtdna_window = is_mtdna_chrom_norm(chrom_norm);
    let is_rdna_contig_window = is_rdna_chrom_norm(chrom_norm, ctx.rdna_contigs);
    let maybe_rdna_intervals = ctx
        .rdna_intervals
        .and_then(|index| index.chroms.get(chrom_norm));

    let region: noodles::core::Region = format!("{}:{}-{}", chrom, window.start + 1, window.end)
        .parse()
        .with_context(|| {
            format!(
                "invalid BAM query region {}:{}-{}",
                chrom,
                window.start + 1,
                window.end
            )
        })?;

    let win_start = window.start as u64;
    let win_end = window.end as u64;
    let mut feature_cursor = maybe_feature.map(|chrom_index| chrom_index.cursor_at(win_start));
    let mut rdna_cursor = maybe_rdna_intervals.map(|chrom_index| chrom_index.cursor_at(win_start));

    let query = reader.query(ctx.header, &region).with_context(|| {
        format!(
            "query failed for region {}:{}-{}",
            chrom,
            window.start + 1,
            window.end
        )
    })?;

    for result in query {
        let record = result.with_context(|| {
            format!(
                "query record error in {}:{}-{}",
                chrom, window.start, window.end
            )
        })?;
        let Some(pos) = alignment_start_0(&record) else {
            continue;
        };
        let owns_start = window_owns_record_start(pos, win_start, win_end);

        if owns_start {
            state.records_seen += 1;
        }
        let flags = record.flags();

        if flags.is_unmapped() {
            if owns_start {
                state.fail_unmapped += 1;
            }
            continue;
        }
        if flags.is_secondary() || flags.is_supplementary() {
            if owns_start {
                state.fail_secondary += 1;
            }
            continue;
        }
        if flags.is_qc_fail() || flags.is_duplicate() {
            if owns_start {
                state.fail_qc += 1;
            }
            continue;
        }
        let record_mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if record_mapq < ctx.config.mapq {
            if owns_start {
                state.fail_mapq += 1;
            }
            continue;
        }

        if ctx.config.r2_only && !flags.is_last_segment() {
            continue;
        }

        let read_match_span = match_span(&record);
        if owns_start {
            state.aligned_qc_reads += 1;
            if is_mtdna_window {
                state.mtdna_reads += 1;
            }

            let mut is_rdna = is_rdna_contig_window;
            if !is_rdna {
                if let Some(chrom_index) = maybe_rdna_intervals {
                    if let Some(cursor) = rdna_cursor.as_mut() {
                        if let Some((first_match, last_match)) = read_match_span {
                            let mid = first_match + (last_match - first_match) / 2;
                            is_rdna = chrom_index.contains(mid, cursor);
                        }
                    }
                }
            }
            if is_rdna {
                state.rdna_reads += 1;
            }
        }

        if let Some(chrom_index) = maybe_feature {
            if let Some((first_match, last_match)) = read_match_span {
                let mid = first_match + (last_match - first_match) / 2;
                if mid >= win_start && mid < win_end {
                    state.total_tags += 1;
                    if let Some(cursor) = feature_cursor.as_mut() {
                        let region = chrom_index.classify(mid, cursor);
                        state.read_dist_counts[region as usize] += 1;
                    }
                }
            }
        }

        if let Some(dense) = maybe_dense {
            if needs_coverage(ctx.config) {
                walk_coverage_record(
                    &record,
                    pos,
                    dense,
                    ctx.coverage_index,
                    ctx.config,
                    state,
                    Some((win_start, win_end)),
                );
            }
        }
    }

    Ok(())
}

fn scan_sequential_bam(ctx: &ScanContext<'_>) -> Result<RnaWorkerState> {
    let file = File::open(ctx.bam_path)?;
    let mut reader = bam::io::Reader::new(file);
    let mut state = RnaWorkerState::new(ctx.config.qc_sample_size);

    let ref_metadata: Vec<_> = ctx
        .header
        .reference_sequences()
        .iter()
        .map(|(name, _)| {
            let chrom = String::from_utf8_lossy(name.as_ref()).to_string();
            let chrom_norm = normalize_chrom(&chrom).into_owned();
            let maybe_dense = ctx.dense_maps.get(&chrom_norm);
            let maybe_feature = ctx
                .feature_index
                .and_then(|idx| idx.chroms.get(&chrom_norm));
            let is_mtdna = is_mtdna_chrom_norm(&chrom_norm);
            let is_rdna_contig = is_rdna_chrom_norm(&chrom_norm, ctx.rdna_contigs);
            let maybe_rdna_intervals = ctx
                .rdna_intervals
                .and_then(|idx| idx.chroms.get(&chrom_norm));
            (
                chrom,
                chrom_norm,
                maybe_dense,
                maybe_feature,
                is_mtdna,
                is_rdna_contig,
                maybe_rdna_intervals,
            )
        })
        .collect();

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {pos} records scanned")?,
    );

    for result in reader.records() {
        let record = result?;
        pb.inc(1);
        state.records_seen += 1;

        let flags = record.flags();
        if flags.is_unmapped() {
            state.fail_unmapped += 1;
            continue;
        }
        if flags.is_secondary() || flags.is_supplementary() {
            state.fail_secondary += 1;
            continue;
        }
        if flags.is_qc_fail() || flags.is_duplicate() {
            state.fail_qc += 1;
            continue;
        }

        let record_mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if record_mapq < ctx.config.mapq {
            state.fail_mapq += 1;
            continue;
        }
        if ctx.config.r2_only && !flags.is_last_segment() {
            continue;
        }

        let id = match record.reference_sequence_id() {
            Some(Ok(id)) => usize::from(id),
            _ => continue,
        };
        let Some((
            _chrom,
            _chrom_norm,
            maybe_dense,
            maybe_feature,
            is_mtdna,
            is_rdna_contig,
            maybe_rdna_intervals,
        )) = ref_metadata.get(id)
        else {
            continue;
        };

        let Some(pos) = alignment_start_0(&record) else {
            continue;
        };

        let read_match_span = match_span(&record);
        state.aligned_qc_reads += 1;
        if *is_mtdna {
            state.mtdna_reads += 1;
        }

        let mut is_rdna = *is_rdna_contig;
        if !is_rdna {
            if let Some(chrom_index) = maybe_rdna_intervals {
                if let Some((first_match, last_match)) = read_match_span {
                    let mid = first_match + (last_match - first_match) / 2;
                    let mut cursor = chrom_index.cursor_at(first_match);
                    is_rdna = chrom_index.contains(mid, &mut cursor);
                }
            }
        }
        if is_rdna {
            state.rdna_reads += 1;
        }

        if let Some(chrom_index) = maybe_feature {
            classify_read_distribution(&mut state, read_match_span, chrom_index, None);
        }

        if let Some(dense) = maybe_dense {
            if needs_coverage(ctx.config) {
                walk_coverage_record(
                    &record,
                    pos,
                    dense,
                    ctx.coverage_index,
                    ctx.config,
                    &mut state,
                    None,
                );
            }
        }
    }
    pb.finish_and_clear();
    Ok(state)
}

fn classify_read_distribution(
    state: &mut RnaWorkerState,
    read_match_span: Option<(u64, u64)>,
    chrom_index: &Arc<ChromFeatureIndex>,
    window: Option<(u64, u64)>,
) {
    if let Some((first_match, last_match)) = read_match_span {
        let mid = first_match + (last_match - first_match) / 2;
        if let Some((win_start, win_end)) = window {
            if mid < win_start || mid >= win_end {
                return;
            }
        }
        let mut cursor = chrom_index.cursor_at(mid);
        let region = chrom_index.classify(mid, &mut cursor);
        state.read_dist_counts[region as usize] += 1;
        state.total_tags += 1;
    }
}

fn needs_coverage(config: &RnaQcConfig) -> bool {
    config.analysis.contains(&AnalysisType::GeneBody)
        || config.analysis.contains(&AnalysisType::ThreePrime)
}

fn walk_coverage_record(
    record: &bam::Record,
    pos: u64,
    dense: &DenseMap,
    index: &AnnotationIndex,
    config: &RnaQcConfig,
    state: &mut RnaWorkerState,
    window: Option<(u64, u64)>,
) {
    if config.ends {
        let is_reverse = record.flags().is_reverse_complemented();
        let target_pos = if is_reverse {
            reference_span(record)
                .map(|(_, end)| end.saturating_sub(1))
                .unwrap_or(pos)
        } else {
            pos
        };
        if let Some((win_start, win_end)) = window {
            if target_pos < win_start || target_pos >= win_end {
                return;
            }
        }
        process_coverage_pos(target_pos, dense, index, config, state);
    } else {
        for_each_aligned_block(record, |block_start, block_end| {
            let (p_lo, p_hi) = if let Some((win_start, win_end)) = window {
                if block_start >= win_end {
                    return;
                }
                let p_lo = block_start.max(win_start);
                let p_hi = block_end.min(win_end);
                if p_lo >= p_hi {
                    return;
                }
                (p_lo, p_hi)
            } else {
                (block_start, block_end)
            };

            if config.step_size > 1 {
                let step = config.step_size as u64;
                let first_sample = if p_lo % step == 0 {
                    p_lo
                } else {
                    p_lo + (step - (p_lo % step))
                };
                for p in (first_sample..p_hi).step_by(config.step_size) {
                    process_coverage_pos(p, dense, index, config, state);
                }
            } else {
                for p in p_lo..p_hi {
                    process_coverage_pos(p, dense, index, config, state);
                }
            }
        });
    }
}

fn process_coverage_pos(
    p: u64,
    dense: &DenseMap,
    index: &AnnotationIndex,
    config: &RnaQcConfig,
    state: &mut RnaWorkerState,
) {
    match dense.get_hits(p) {
        Hits::None => {}
        Hits::Single(g_idx) => process_gene_hit(g_idx as usize, p, index, config, state),
        Hits::Multi(indices) => {
            for &g_idx in indices {
                process_gene_hit(g_idx as usize, p, index, config, state);
            }
        }
    }
}

fn process_gene_hit(
    g_idx: usize,
    p: u64,
    index: &AnnotationIndex,
    config: &RnaQcConfig,
    state: &mut RnaWorkerState,
) {
    state.overlaps_found += 1;
    let Some(gene) = index.genes.get(g_idx) else {
        return;
    };
    if let Some(s_5p) = gene.bin_map.get_spliced_5p(p) {
        let pct_idx = (s_5p as usize * 100) / (gene.total_len as usize).max(1);
        gene.add_percentile(pct_idx, 1);

        let dist_3p = (gene.total_len as usize).saturating_sub(s_5p as usize + 1);
        let bin = dist_3p / config.three_prime_bin_size;
        gene.add_3p(bin, 1);
    }
}
