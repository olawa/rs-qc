mod aggregate;
mod config;
mod report;
mod scan;
pub(crate) mod state;

use crate::analysis::alignment_qc::sample_name_from_alignment_path;
use crate::analysis::contamination::ContaminantIndex;
use crate::analysis::feature_index::FeatureIndex;
use crate::analysis::index::{AnnotationIndex, DenseMap};
use crate::analysis::qc::{InlineQcState, ReadEndType};
use crate::analysis::types::AnalysisType;
use crate::io::annotation::{load_annotation, load_genes, AnnotationConfig, AnnotationFormat};
use crate::io::bam::alignment_start_0;
use crate::models::{normalize_chrom, Gene};
use crate::stats::plotting::PlotMetadata;
use anyhow::{bail, Context, Result};
use noodles::bam;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::num::NonZeroUsize;
use std::sync::Arc;
use indicatif::ProgressBar;
use rayon::prelude::*;
use std::time::Instant;

pub use config::{DenseMapScope, RnaQcConfig};

pub fn run_rna(config: RnaQcConfig) -> Result<()> {
    let threads = normalize_thread_count(config.threads);
    if config.three_prime_bin_size == 0 {
        bail!("--three-prime-bin-size must be greater than 0");
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap_or_else(|_| {}); // May already be initialized

    let start_time = Instant::now();

    let ann_config = AnnotationConfig {
        gene_id_delimiter: config.gene_id_delimiter,
        gene_id_regex: config.gene_id_regex.clone(),
        biotype_filter: Some(config.biotype.clone()),
        three_prime_cluster_window: config.three_prime_cluster_window,
        min_transcript_length: config.min_length,
        max_3p_dist: config.max_3p_dist,
        transcript_centric: config.transcript_centric,
        plus_strand_only: config.plus,
        isoform_select: config.isoform_select.clone(),
        strict_cluster: config.strict_cluster,
        three_prime_bin_size: config.three_prime_bin_size,
    };
    let ann_format = match config.annotation_format.to_lowercase().as_str() {
        "bed12" => AnnotationFormat::Bed12,
        "gtf" => AnnotationFormat::Gtf,
        _ => AnnotationFormat::Auto,
    };

    println!("--------------------------------------------------");
    println!("rs-qc rna v{}", env!("CARGO_PKG_VERSION"));
    println!("--------------------------------------------------");
    println!("Loading annotation: {}", config.annotation);

    let index_path = format!("{}.ridx", config.annotation);
    let dist_enabled = config.analysis.contains(&AnalysisType::Distribution);
    let (index_proto, feature_index) = if dist_enabled {
        if config.load_index {
            println!(
                "  - Distribution analysis needs transcript annotations, so the source annotation will be parsed instead of loading only the collapsed index."
            );
        }
        println!("  - Parsing source annotation: {}", config.annotation);
        let loaded = load_annotation(&config.annotation, ann_format, &ann_config)?;
        println!("  - Extracted {} eligible genes.", loaded.genes.len());
        println!("  - Building distribution feature index...");
        let feature_index = Arc::new(FeatureIndex::build(&loaded.transcripts));
        println!("  - Initializing annotation index...");
        let index = AnnotationIndex::new(loaded.genes, false);
        if config.save_index {
            println!("  - Saving index for future use: {}", index_path);
            index.save_to_file(&index_path)?;
        }
        (index, Some(feature_index))
    } else if config.load_index && std::path::Path::new(&index_path).exists() {
        println!("  - Loading pre-built index: {}", index_path);
        (
            AnnotationIndex::load_from_file(
                &index_path,
                config.max_3p_dist,
                config.three_prime_bin_size,
            )?,
            None,
        )
    } else {
        let genes = load_genes(&config.annotation, ann_format, &ann_config)?;
        println!("  - Extracted {} eligible genes.", genes.len());
        let index = AnnotationIndex::new(genes, false);

        if config.save_index {
            println!("  - Saving index for future use: {}", index_path);
            index.save_to_file(&index_path)?;
        }
        (index, None)
    };

    let needs_coverage = config.analysis.contains(&AnalysisType::GeneBody)
        || config.analysis.contains(&AnalysisType::ThreePrime)
        || config.analysis.contains(&AnalysisType::Qc);
    let coverage_index = if needs_coverage {
        Some(autosomal_coverage_index(&index_proto))
    } else {
        None
    };
    let coverage_index_ref = coverage_index.as_ref().unwrap_or(&index_proto);

    println!(
        "  - Index ready with {} active genes.",
        coverage_index_ref.genes.len()
    );

    let rdna_contigs = Arc::new(normalized_name_set(&config.rdna_contigs));
    let rdna_intervals = if let Some(path) = &config.rdna_bed {
        println!("  - Loading rDNA intervals: {}", path);
        Some(Arc::new(ContaminantIndex::from_bed(path)?))
    } else {
        None
    };

    let mut classic_all: HashMap<String, Vec<f64>> = HashMap::new();
    let mut classic_percent_all: HashMap<String, Vec<f64>> = HashMap::new();
    let mut dist_3p_all: HashMap<String, Vec<f64>> = HashMap::new();
    let mut dist_3p_support_all: HashMap<String, Vec<usize>> = HashMap::new();
    let mut dist_3p_raw_all: HashMap<String, Vec<f64>> = HashMap::new();
    let mut dist_3p_raw_analyzed_all: HashMap<String, Vec<f64>> = HashMap::new();
    let mut inner_distance_all: HashMap<String, Vec<f64>> = HashMap::new();
    let mut plot_metadata_all: HashMap<String, PlotMetadata> = HashMap::new();

    for (i, bam_path) in config.input.iter().enumerate() {
        let sample_name = sample_name_from_alignment_path(bam_path);
        println!(
            "\n[{}/{}] Processing BAM: {}",
            i + 1,
            config.input.len(),
            bam_path
        );

        let header = {
            let file = File::open(bam_path)?;
            let mut reader = bam::io::Reader::new(file);
            reader.read_header()?
        };

        let bai_path = crate::analysis::bam_scan::find_bai_path(bam_path);
        let mut total_state = if let (Some(bp), scope) = (&bai_path, config.dense_map_scope) {
            let scan_start = Instant::now();
            let res = match scope {
                DenseMapScope::Window => {
                    println!(
                        "  - DenseMap scope: Window, window size: {} bp",
                        config.window_size
                    );
                    let windows =
                        crate::analysis::bam_scan::generate_windows(&header, config.window_size);
                    let pb = ProgressBar::new(windows.len() as u64);
                    pb.set_style(indicatif::ProgressStyle::default_bar().template(
                        "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Windows ({eta})",
                    )?);

                    let res = scan::scan_windows_with_pb(
                        &scan::ScanContext {
                            bam_path,
                            header: &header,
                            dense_maps: &HashMap::new(),
                            coverage_index: coverage_index_ref,
                            feature_index: feature_index.as_deref(),
                            rdna_contigs: rdna_contigs.as_ref(),
                            rdna_intervals: rdna_intervals.as_deref(),
                            config: &config,
                            threads,
                        },
                        bp,
                        &windows,
                        Some(pb.clone()),
                    )?;
                    pb.finish_and_clear();
                    res
                }
                DenseMapScope::Chunk | DenseMapScope::Chrom => {
                    let mut state = state::RnaWorkerState::new(&config, config.qc_sample_size);
                    let chunk_size = if config.dense_map_scope == DenseMapScope::Chrom {
                        u64::MAX
                    } else {
                        config.dense_map_chunk_size
                    };

                    println!(
                        "  - DenseMap scope: {:?}, chunk size: {} bp",
                        config.dense_map_scope, chunk_size
                    );

                    let all_windows =
                        crate::analysis::bam_scan::generate_windows(&header, config.window_size);
                    let pb = ProgressBar::new(all_windows.len() as u64);
                    pb.set_style(indicatif::ProgressStyle::default_bar().template(
                        "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Windows ({eta})",
                    )?);

                    for (name, seq) in header.reference_sequences() {
                        let chrom = String::from_utf8_lossy(name.as_ref()).to_string();
                        let chrom_norm = normalize_chrom(&chrom).into_owned();
                        let chrom_len = seq.length().get() as u64;

                        let (start, end) = coverage_index_ref
                            .chrom_spans
                            .get(&chrom_norm)
                            .cloned()
                            .unwrap_or((0, chrom_len));

                        let mut chunk_start = start;
                        while chunk_start < end {
                            let chunk_end = (chunk_start + chunk_size).min(end);
                            let dense = coverage_index_ref.build_dense_map_for_range(
                                &chrom,
                                chunk_start,
                                chunk_end,
                            );

                            let mut dense_maps = HashMap::new();
                            if let Some(d) = dense {
                                dense_maps.insert(chrom_norm.clone(), Arc::new(d));
                            }

                            let windows = crate::analysis::bam_scan::generate_windows_for_range(
                                &chrom,
                                &chrom_norm,
                                chunk_start,
                                chunk_end,
                                config.window_size,
                            );

                            let chunk_state = scan::scan_windows_with_pb(
                                &scan::ScanContext {
                                    bam_path,
                                    header: &header,
                                    dense_maps: &dense_maps,
                                    coverage_index: coverage_index_ref,
                                    feature_index: feature_index.as_deref(),
                                    rdna_contigs: rdna_contigs.as_ref(),
                                    rdna_intervals: rdna_intervals.as_deref(),
                                    config: &config,
                                    threads,
                                },
                                bp,
                                &windows,
                                Some(pb.clone()),
                            )?;

                            state = state.merge(chunk_state);
                            chunk_start = chunk_end;
                        }
                    }
                    pb.finish_and_clear();
                    state
                }
                DenseMapScope::All => {
                    let dense_maps = build_all_dense_maps(coverage_index_ref, &config)?;
                    scan::scan_bam_coverage(&scan::ScanContext {
                        bam_path,
                        header: &header,
                        dense_maps: &dense_maps,
                        coverage_index: coverage_index_ref,
                        feature_index: feature_index.as_deref(),
                        rdna_contigs: rdna_contigs.as_ref(),
                        rdna_intervals: rdna_intervals.as_deref(),
                        config: &config,
                        threads,
                    })?
                }
            };
            println!("  - Main BAM scan took: {:?}", scan_start.elapsed());
            res
        } else {
            // Fallback to All mode if no BAI or explicitly requested
            if bai_path.is_none() && config.dense_map_scope != DenseMapScope::All {
                println!("  - Warning: No BAI index found. Falling back to 'All' DenseMap scope for sequential scan.");
            }

            let dense_maps = build_all_dense_maps(coverage_index_ref, &config)?;
            scan::scan_bam_coverage(&scan::ScanContext {
                bam_path,
                header: &header,
                dense_maps: &dense_maps,
                coverage_index: coverage_index_ref,
                feature_index: feature_index.as_deref(),
                rdna_contigs: rdna_contigs.as_ref(),
                rdna_intervals: rdna_intervals.as_deref(),
                config: &config,
                threads,
            })?
        };

        if config.analysis.contains(&AnalysisType::Qc) {
            let qc_start = Instant::now();
            total_state.qc = scan_inline_qc_sample(
                bam_path,
                &header,
                coverage_index_ref,
                &config,
            )?;
            println!("  - Read-name pair QC scan took: {:?}", qc_start.elapsed());
        }

        let finalize_start = Instant::now();
        let passed_filters = total_state.records_seen
            - total_state.fail_unmapped
            - total_state.fail_secondary
            - total_state.fail_qc
            - total_state.fail_mapq;

        println!(
            "  - Results: {} exonic overlaps processed from {} records.",
            total_state.overlaps_found, total_state.records_seen
        );
        println!("    - Passed filters:   {}", passed_filters);
        println!(
            "    - Failed MapQ < {}: {}",
            config.mapq, total_state.fail_mapq
        );

        let stats_res = Arc::new(std::sync::Mutex::new(None));

        let report_errors = Arc::new(std::sync::Mutex::new(Vec::<String>::new()));

        rayon::scope(|s| {
            if config.analysis.contains(&AnalysisType::GeneBody)
                || config.analysis.contains(&AnalysisType::ThreePrime)
            {
                let stats_ptr = Arc::clone(&stats_res);
                let errors = Arc::clone(&report_errors);
                let config_ref = &config;
                let sample_ref = &sample_name;
                let state_ref = &total_state;
                let index_ref = coverage_index_ref;
                s.spawn(move |_| {
                    let agg_start = Instant::now();
                    let stats = aggregate::aggregate_sample(
                        index_ref,
                        config_ref,
                        state_ref,
                    );
                    println!("  - Coverage aggregation took: {:?}", agg_start.elapsed());

                    if let Err(e) =
                        report::write_sample_gene_body_plot(config_ref, sample_ref, &stats)
                    {
                        errors.lock().unwrap().push(format!(
                            "failed to write gene body plot for {sample_ref}: {e}"
                        ));
                    }
                    *stats_ptr.lock().unwrap() = Some(stats);
                });
            }

            if config.analysis.contains(&AnalysisType::Qc) {
                let errors = Arc::clone(&report_errors);
                let config_ref = &config;
                let sample_ref = &sample_name;
                let state_ref = &total_state;
                s.spawn(move |_| {
                    let agg_start = Instant::now();
                    let mut qc = state_ref.qc.summary.clone();
                    qc.aligned_qc_reads = state_ref.aligned_qc_reads;
                    qc.mtdna_reads = state_ref.mtdna_reads;
                    qc.rdna_reads = state_ref.rdna_reads;

                    let qc_path = format!("{}.{}.rna_qc.txt", config_ref.output, sample_ref);
                    if let Err(e) = qc.write_summary_file(&qc_path) {
                        errors.lock().unwrap().push(format!(
                            "failed to write RNA QC summary for {sample_ref}: {e}"
                        ));
                    }

                    println!("  - QC metrics aggregation took: {:?}", agg_start.elapsed());
                    println!("{}", qc.summary_text());
                });
            }

            if config.analysis.contains(&AnalysisType::Distribution) {
                let errors = Arc::clone(&report_errors);
                let config_ref = &config;
                let sample_ref = &sample_name;
                let state_ref = &total_state;
                let f_idx_clone = feature_index.clone();
                s.spawn(move |_| {
                    let agg_start = Instant::now();
                    if let Some(f_idx) = f_idx_clone.as_ref() {
                        let dist_path =
                            format!("{}.{}.read_distribution.txt", config_ref.output, sample_ref);
                        if let Err(e) = state_ref.write_read_distribution_report(&dist_path, f_idx)
                        {
                            errors.lock().unwrap().push(format!(
                                "failed to write read distribution for {sample_ref}: {e}"
                            ));
                        }
                    }
                    println!("  - Distribution analysis took: {:?}", agg_start.elapsed());
                });
            }
        });

        let errors = report_errors.lock().unwrap();
        if !errors.is_empty() {
            bail!("{}", errors.join("; "));
        }

        // Update multi-sample maps safely outside the scope
        let stats_opt = stats_res.lock().unwrap().take();
        if let Some(stats) = stats_opt {
            let mut visual_qc = total_state.qc.summary.clone();
            visual_qc.aligned_qc_reads = total_state.aligned_qc_reads;
            visual_qc.mtdna_reads = total_state.mtdna_reads;
            visual_qc.rdna_reads = total_state.rdna_reads;

            let report_start = Instant::now();
            report::write_sample_reports(&config, &sample_name, &stats)?;
            report::write_sample_visual_reports(
                &config,
                bam_path,
                &sample_name,
                &stats,
                &visual_qc,
                &total_state.read_dist_counts,
                index_proto.genes.as_slice(),
                ann_format,
            )?;
            println!("  - Sample-specific reports and plots took: {:?}", report_start.elapsed());
            classic_all.insert(sample_name.to_string(), stats.percentile_means.clone());
            classic_percent_all.insert(sample_name.to_string(), stats.percentile_normalized);
            dist_3p_all.insert(sample_name.to_string(), stats.dist_3p_means.clone());
            dist_3p_support_all.insert(sample_name.to_string(), stats.dist_3p_support.clone());
            dist_3p_raw_all.insert(sample_name.to_string(), stats.dist_3p_sums_raw_all.clone());
            dist_3p_raw_analyzed_all.insert(
                sample_name.to_string(),
                stats.dist_3p_sums_raw_analyzed.clone(),
            );
            plot_metadata_all.insert(
                sample_name.to_string(),
                PlotMetadata {
                    total_reads: stats.total_reads,
                    active_genes: stats.active_genes,
                    active_3p_genes: stats.active_3p_genes,
                },
            );
        }

        if config.analysis.contains(&AnalysisType::Qc) {
            let mut qc = total_state.qc.summary.clone();
            qc.aligned_qc_reads = total_state.aligned_qc_reads;
            qc.mtdna_reads = total_state.mtdna_reads;
            qc.rdna_reads = total_state.rdna_reads;
            let series: Vec<f64> = qc
                .clipped_inner_distance_series(-200, 200)
                .into_iter()
                .map(|v| v as f64)
                .collect();
            inner_distance_all.insert(sample_name.to_string(), series);
        }

        println!(
            "  - Total aggregation and reporting took: {:?}",
            finalize_start.elapsed()
        );
        coverage_index_ref.reset_coverage();
    }

    let multi_start = Instant::now();
    report::write_multi_sample_outputs(
        &config,
        &classic_all,
        &classic_percent_all,
        &dist_3p_all,
        &dist_3p_support_all,
        &dist_3p_raw_all,
        &dist_3p_raw_analyzed_all,
        &inner_distance_all,
        &plot_metadata_all,
    )?;
    println!("  - Multi-sample aggregation and TSV/SVG generation took: {:?}", multi_start.elapsed());

    println!("\\nAll tasks completed in {:.2?}.", start_time.elapsed());
    Ok(())
}

fn normalized_name_set(names: &[String]) -> HashSet<String> {
    scan::normalized_name_set(names)
}

fn autosomal_coverage_index(index: &AnnotationIndex) -> AnnotationIndex {
    scan::autosomal_coverage_index(index)
}

fn normalize_thread_count(threads: usize) -> usize {
    scan::normalize_thread_count(threads)
}

fn build_all_dense_maps(
    index: &AnnotationIndex,
    config: &RnaQcConfig,
) -> Result<HashMap<String, Arc<DenseMap>>> {
    let memory_estimate_bytes = index.estimate_dense_map_memory();
    println!(
        "  - Pre-calculating all dense mappings (est. peak memory: {:.2} GB)...",
        memory_estimate_bytes as f64 / 1_073_741_824.0
    );

    let index_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.dense_map_workers.max(1))
        .build()
        .with_context(|| "failed to create index thread pool")?;

    let dense_maps: HashMap<String, Arc<DenseMap>> = index_pool.install(|| {
        index
            .chrom_spans
            .keys()
            .collect::<Vec<_>>()
            .into_par_iter()
            .filter_map(|chrom| {
                index
                    .build_dense_map(chrom)
                    .map(|dense| (chrom.to_string(), Arc::new(dense)))
            })
            .fold(
                || HashMap::new(),
                |mut d_acc, (chrom, dense)| {
                    d_acc.insert(chrom, dense);
                    d_acc
                },
            )
            .reduce(
                || HashMap::new(),
                |mut d1, d2| {
                    d1.extend(d2);
                    d1
                },
            )
    });

    Ok(dense_maps)
}

#[cfg(test)]
fn window_owns_record_start(pos: u64, win_start: u64, win_end: u64) -> bool {
    scan::window_owns_record_start(pos, win_start, win_end)
}

fn scan_inline_qc_sample(
    bam_path: &str,
    header: &noodles::sam::Header,
    index: &AnnotationIndex,
    config: &RnaQcConfig,
) -> Result<InlineQcState> {
    let file = File::open(bam_path)?;
    let worker_count = NonZeroUsize::new(normalize_thread_count(config.threads)).unwrap();
    let mt_reader = noodles::bgzf::MultithreadedReader::with_worker_count(worker_count, file);
    let mut reader = bam::io::Reader::from(mt_reader);
    let _ = reader.read_header()?;
    let mut qc = InlineQcState::new(config.qc_sample_size);

    let ref_metadata: Vec<_> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| {
            let chrom = String::from_utf8_lossy(name.as_ref()).to_string();
            let chrom_norm = normalize_chrom(&chrom).into_owned();
            (chrom, chrom_norm)
        })
        .collect();

    let mut current_dense: Option<(String, Arc<DenseMap>)> = None;

    for result in reader.records() {
        if !qc.needs_more() {
            break;
        }
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_qc_fail()
            || flags.is_duplicate()
        {
            continue;
        }

        let record_mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if record_mapq < config.mapq {
            continue;
        }
        if config.r2_only && !flags.is_last_segment() {
            continue;
        }

        let id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            _ => continue,
        };
        let Some((chrom, chrom_norm)) = ref_metadata.get(id) else {
            continue;
        };

        // On-demand DenseMap building
        if current_dense
            .as_ref()
            .map(|(c, _)| c != chrom_norm)
            .unwrap_or(true)
        {
            if let Some(dense) = index.build_dense_map(chrom) {
                current_dense = Some((chrom_norm.clone(), Arc::new(dense)));
            } else {
                continue;
            }
        }

        let dense = &current_dense.as_ref().unwrap().1;
        maybe_observe_inline_qc(&mut qc, &record, chrom, dense, &index.genes, 0, u64::MAX);
    }

    Ok(qc)
}

fn maybe_observe_inline_qc(
    qc: &mut InlineQcState,
    record: &bam::Record,
    chrom: &str,
    dense: &DenseMap,
    genes: &[Gene],
    win_start: u64,
    win_end: u64,
) {
    let Some(pos) = alignment_start_0(record) else {
        return;
    };
    if pos >= win_start && pos < win_end {
        let is_reverse = record.flags().is_reverse_complemented();
        let end_type = if is_reverse {
            ReadEndType::ThreePrime
        } else {
            ReadEndType::FivePrime
        };
        qc.observe_read_end(record, chrom, dense, genes, end_type);

        // Inner distance is already observed inside observe_read_end
    }
}

pub fn default_rdna_contigs() -> Vec<String> {
    vec![
        "rdna".to_string(),
        "rrna".to_string(),
        "45s".to_string(),
        "18s".to_string(),
        "28s".to_string(),
        "5s".to_string(),
        "rn45s".to_string(),
        "rn18s".to_string(),
        "rn28s".to_string(),
        "rn5s".to_string(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_threads_coerces_zero_to_one() {
        assert_eq!(normalize_thread_count(0), 1);
        assert_eq!(normalize_thread_count(4), 4);
    }

    #[test]
    fn window_owns_record_start_only_inside_half_open_window() {
        assert!(!window_owns_record_start(9, 10, 20));
        assert!(window_owns_record_start(10, 10, 20));
        assert!(window_owns_record_start(19, 10, 20));
        assert!(!window_owns_record_start(20, 10, 20));
    }
}
