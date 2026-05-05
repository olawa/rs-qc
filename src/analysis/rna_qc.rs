mod aggregate;
mod config;
mod report;
mod scan;
mod state;

use crate::analysis::alignment_qc::sample_name_from_alignment_path;
use crate::analysis::bam_scan::{find_bai_path, generate_windows};
use crate::analysis::contamination::ContaminantIndex;
use crate::analysis::feature_index::FeatureIndex;
use crate::analysis::index::{AnnotationIndex, DenseMap};
use crate::analysis::qc::{InlineQcState, ReadEndType};
use crate::analysis::types::AnalysisType;
use crate::io::annotation::{load_annotation, load_genes, AnnotationConfig, AnnotationFormat};
use crate::io::bam::{alignment_start_0, for_each_aligned_block, match_span, reference_span};
use crate::models::{normalize_chrom, Gene};
use crate::stats::plotting::PlotMetadata;
use anyhow::{bail, Result};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use noodles::bam;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::num::NonZeroUsize;
use std::sync::Arc;
use std::time::Instant;

thread_local! {
    static BAM_READER: RefCell<Option<(String, bam::io::IndexedReader<noodles::bgzf::Reader<File>>)>> = RefCell::new(None);
}

const MAX_REASONABLE_INNER_DISTANCE_BP: i64 = 2_000;
pub use config::RnaQcConfig;
use state::RnaWorkerState;

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

        println!(
            "  - Pre-calculating metadata and dense mappings for {} chromosomes...",
            coverage_index_ref.chrom_spans.len()
        );
        use rayon::prelude::*;
        let dense_maps: HashMap<String, Arc<DenseMap>> = coverage_index_ref
            .chrom_spans
            .keys()
            .collect::<Vec<_>>()
            .into_par_iter()
            .filter_map(|chrom| {
                coverage_index_ref
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
            );

        let dense_maps = Arc::new(dense_maps);

        let mb_index_path = find_bai_path(bam_path);
        let finalize_start = Instant::now();
        let mut total_state = if let Some(bai_path) = mb_index_path {
            println!("  - BAI Index found. Using high-performance parallel dense scan.");
            let bai = bam::bai::read(&bai_path)?;
            let windows = generate_windows(&header, 10_000_000);

            let m_pb = MultiProgress::new();
            let pb = m_pb.add(ProgressBar::new(windows.len() as u64));
            pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Windows ({eta})")?);

            let bai_arc = Arc::new(bai);
            let n_genes = coverage_index_ref.genes.len();
            let max_3p_dist = config.max_3p_dist;
            let qc_sample_size = config.qc_sample_size.div_ceil(threads);
            let total_state = windows
                .par_iter()
                .fold(
                    move || RnaWorkerState::new_with_capacity(n_genes, max_3p_dist, qc_sample_size),
                    |mut state, window| {
                    BAM_READER.with(|cell| {
                        let mut opt = cell.borrow_mut();
                        let needs_reader = opt
                            .as_ref()
                            .map(|(path, _)| path != bam_path)
                            .unwrap_or(true);
                        if needs_reader {
                            let f = File::open(bam_path).expect("Failed to open BAM");
                            let mut r = bam::io::indexed_reader::Builder::default()
                                .set_index(bai_arc.as_ref().clone())
                                .build_from_reader(f)
                                .expect("Failed to build indexed reader");
                            let _ = r.read_header().expect("Failed to read BAM header");
                            *opt = Some((bam_path.clone(), r));
                        }

                        let reader = &mut opt.as_mut().unwrap().1;
                        let chrom = &window.chrom;
                        let chrom_norm = &window.chrom_norm;
                        let maybe_dense = dense_maps.get(chrom_norm);
                        let maybe_feature = feature_index
                            .as_ref()
                            .and_then(|index| index.chroms.get(chrom_norm));
                        let is_mtdna_window = is_mtdna_chrom_norm(chrom_norm);
                        let is_rdna_contig_window = is_rdna_chrom_norm(chrom_norm, rdna_contigs.as_ref());
                        let maybe_rdna_intervals = rdna_intervals
                            .as_ref()
                            .and_then(|index| index.chroms.get(chrom_norm));

                        let region: noodles::core::Region =
                            match format!("{}:{}-{}", chrom, window.start + 1, window.end).parse() {
                                Ok(r) => r,
                                Err(_) => return,
                            };

                        let win_start = window.start as u64;
                        let win_end = window.end as u64;
                        let mut feature_cursor = maybe_feature.map(|chrom_index| chrom_index.cursor_at(win_start));
                        let mut rdna_cursor = maybe_rdna_intervals
                            .map(|chrom_index| chrom_index.cursor_at(win_start));

                        match reader.query(&header, &region) {
                            Ok(query) => {
                                for result in query {
                                    match result {
                                        Ok(record) => {
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
                                            let record_mapq = record
                                                .mapping_quality()
                                                .map(|m| m.get())
                                                .unwrap_or(255);
                                            if record_mapq < config.mapq {
                                                if owns_start {
                                                    state.fail_mapq += 1;
                                                }
                                                continue;
                                            }

                                            if config.r2_only && !flags.is_last_segment() {
                                                continue;
                                            }

                                            let genes_ref = &coverage_index_ref.genes;
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
                                                                let mid = first_match
                                                                    + (last_match - first_match) / 2;
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
                                                let needs_coverage =
                                                    config.analysis.contains(&AnalysisType::GeneBody)
                                                        || config
                                                            .analysis
                                                            .contains(&AnalysisType::ThreePrime);

                                                if needs_coverage {
                                                    let process_pos = |p: u64, state: &mut RnaWorkerState| {
                                                        use crate::analysis::index::Hits;
                                                        match dense.get_hits(p) {
                                                            Hits::None => {}
                                                            Hits::Single(g_idx) => {
                                                                let g_idx = g_idx as usize;
                                                                state.overlaps_found += 1;
                                                                let gene = &genes_ref[g_idx];
                                                                if let Some(s_5p) = gene.bin_map.get_spliced_5p(p) {
                                                                    let pct_idx = (s_5p as usize * 100) / (gene.total_len as usize).max(1);
                                                                    gene.add_percentile(pct_idx, 1);

                                                                    let dist_3p = (gene.total_len as usize).saturating_sub(s_5p as usize + 1);
                                                                    let bin = dist_3p / config.three_prime_bin_size;
                                                                    gene.add_3p(bin, 1);
                                                                }
                                                            }
                                                            Hits::Multi(indices) => {
                                                                for &g_idx in indices {
                                                                    let g_idx = g_idx as usize;
                                                                    state.overlaps_found += 1;
                                                                        let gene = &genes_ref[g_idx];
                                                                    if let Some(s_5p) = gene.bin_map.get_spliced_5p(p) {
                                                                        let pct_idx = (s_5p as usize * 100) / (gene.total_len as usize).max(1);
                                                                        gene.add_percentile(pct_idx, 1);

                                                                        let dist_3p = (gene.total_len as usize).saturating_sub(s_5p as usize + 1);
                                                                        let bin = dist_3p / config.three_prime_bin_size;
                                                                        gene.add_3p(bin, 1);
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    };

                                                if config.ends {
                                                    let is_reverse = flags.is_reverse_complemented();
                                                    let target_pos = if is_reverse {
                                                        reference_span(&record)
                                                            .map(|(_, end)| end.saturating_sub(1))
                                                            .unwrap_or(pos)
                                                    } else {
                                                        pos
                                                    };
                                                    process_pos(target_pos, &mut state);
                                                } else {
                                                    for_each_aligned_block(&record, |block_start, block_end| {
                                                        if block_start >= win_end {
                                                            return;
                                                        }
                                                        let p_lo = block_start.max(win_start);
                                                        let p_hi = block_end.min(win_end);
                                                        if p_lo >= p_hi {
                                                            return;
                                                        }

                                                        if config.step_size > 1 {
                                                            let step = config.step_size as u64;
                                                            let first_sample = if p_lo % step == 0 {
                                                                p_lo
                                                            } else {
                                                                p_lo + (step - (p_lo % step))
                                                            };
                                                            for p in (first_sample..p_hi).step_by(config.step_size) {
                                                                process_pos(p, &mut state);
                                                            }
                                                        } else {
                                                            for p in p_lo..p_hi {
                                                                process_pos(p, &mut state);
                                                            }
                                                        }
                                                    });
                                                }
                                            }
                                        }
                                    }
                                    Err(e) => {
                                            eprintln!(
                                                "    [ERROR] Query record error in {}:{}-{}: {:?}",
                                                chrom, window.start, window.end, e
                                            );
                                            return;
                                        }
                                    }
                                }
                            }
                            Err(e) => {
                                println!(
                                    "    [ERROR] Query failed for region {}:{}-{}: {}",
                                    chrom,
                                    window.start + 1,
                                    window.end,
                                    e
                                );
                            }
                        }
                    });
                    pb.inc(1);
                    state
                })
                .reduce(
                    move || RnaWorkerState::new_with_capacity(n_genes, max_3p_dist, qc_sample_size),
                    RnaWorkerState::merge,
                );
            total_state
        } else {
            println!("  - No BAI Index found. Falling back to sequential single-threaded scan.");
            let file = File::open(bam_path)?;
            let mut reader = bam::io::Reader::new(file);
            let n_genes = coverage_index_ref.genes.len();
            let mut state = RnaWorkerState::new_with_capacity(
                n_genes,
                config.max_3p_dist,
                config.qc_sample_size,
            );

            let ref_metadata: Vec<_> = header
                .reference_sequences()
                .iter()
                .map(|(name, _)| {
                    let chrom = String::from_utf8_lossy(name.as_ref()).to_string();
                    let chrom_norm = normalize_chrom(&chrom).into_owned();
                    let maybe_dense = dense_maps.get(&chrom_norm);
                    let maybe_feature = feature_index
                        .as_ref()
                        .and_then(|idx| idx.chroms.get(&chrom_norm));
                    let is_mtdna = is_mtdna_chrom_norm(&chrom_norm);
                    let is_rdna_contig = is_rdna_chrom_norm(&chrom_norm, rdna_contigs.as_ref());
                    let maybe_rdna_intervals = rdna_intervals
                        .as_ref()
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
                if record_mapq < config.mapq {
                    state.fail_mapq += 1;
                    continue;
                }
                if config.r2_only && !flags.is_last_segment() {
                    continue;
                }

                let id = match record.reference_sequence_id() {
                    Some(Ok(id)) => usize::from(id),
                    _ => continue,
                };
                let (
                    _chrom,
                    _chrom_norm,
                    maybe_dense,
                    maybe_feature,
                    is_mtdna,
                    is_rdna_contig,
                    maybe_rdna_intervals,
                ) = &ref_metadata[id];

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
                    if let Some((first_match, last_match)) = read_match_span {
                        let mid = first_match + (last_match - first_match) / 2;
                        let mut cursor = chrom_index.cursor_at(mid);
                        let region = chrom_index.classify(mid, &mut cursor);
                        state.read_dist_counts[region as usize] += 1;
                        state.total_tags += 1;
                    }
                }

                if let Some(dense) = maybe_dense {
                    if config.analysis.contains(&AnalysisType::GeneBody)
                        || config.analysis.contains(&AnalysisType::ThreePrime)
                    {
                        let process_pos = |p: u64, state: &mut RnaWorkerState| {
                            use crate::analysis::index::Hits;
                            match dense.get_hits(p) {
                                Hits::None => {}
                                Hits::Single(g_idx) => {
                                    let g_idx = g_idx as usize;
                                    state.overlaps_found += 1;
                                    let gene = &coverage_index_ref.genes[g_idx];
                                    if let Some(s_5p) = gene.bin_map.get_spliced_5p(p) {
                                        let pct_idx = (s_5p as usize * 100)
                                            / (gene.total_len as usize).max(1);
                                        gene.add_percentile(pct_idx, 1);

                                        let dist_3p = (gene.total_len as usize)
                                            .saturating_sub(s_5p as usize + 1);
                                        let bin = dist_3p / config.three_prime_bin_size;
                                        gene.add_3p(bin, 1);
                                    }
                                }
                                Hits::Multi(indices) => {
                                    for &g_idx in indices {
                                        let g_idx = g_idx as usize;
                                        state.overlaps_found += 1;
                                        let gene = &coverage_index_ref.genes[g_idx];
                                        if let Some(s_5p) = gene.bin_map.get_spliced_5p(p) {
                                            let pct_idx = (s_5p as usize * 100)
                                                / (gene.total_len as usize).max(1);
                                            gene.add_percentile(pct_idx, 1);

                                            let dist_3p = (gene.total_len as usize)
                                                .saturating_sub(s_5p as usize + 1);
                                            let bin = dist_3p / config.three_prime_bin_size;
                                            gene.add_3p(bin, 1);
                                        }
                                    }
                                }
                            }
                        };

                        if config.ends {
                            let is_reverse = flags.is_reverse_complemented();
                            let target_pos = if is_reverse {
                                reference_span(&record)
                                    .map(|(_, end)| end.saturating_sub(1))
                                    .unwrap_or(pos)
                            } else {
                                pos
                            };
                            process_pos(target_pos, &mut state);
                        } else {
                            for_each_aligned_block(&record, |block_start, block_end| {
                                if config.step_size > 1 {
                                    for p in (block_start..block_end).step_by(config.step_size) {
                                        process_pos(p, &mut state);
                                    }
                                } else {
                                    for p in block_start..block_end {
                                        process_pos(p, &mut state);
                                    }
                                }
                            });
                        }
                    }
                }
            }
            pb.finish_and_clear();
            state
        };

        if config.analysis.contains(&AnalysisType::Qc) {
            let qc_start = Instant::now();
            total_state.qc = scan_inline_qc_sample(
                bam_path,
                &header,
                dense_maps.as_ref(),
                coverage_index_ref,
                &config,
            )?;
            println!("  - Read-name pair QC scan took: {:?}", qc_start.elapsed());
        }

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
                let _errors = Arc::clone(&report_errors);
                let config_ref = &config;
                let sample_ref = &sample_name;
                let state_ref = &total_state;
                let index_ref = coverage_index_ref;
                s.spawn(move |_| {
                    let agg_start = Instant::now();
                    let stats = aggregate::aggregate_sample(
                        index_ref,
                        config_ref,
                        state_ref.aligned_qc_reads,
                    );
                    println!("  - Coverage aggregation took: {:?}", agg_start.elapsed());

                    let _ = report::write_sample_gene_body_plot(config_ref, sample_ref, &stats);
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
            let _ = report::write_sample_reports(&config, &sample_name, &stats);
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
            let qc = total_state.qc.summary.clone();
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

fn window_owns_record_start(pos: u64, win_start: u64, win_end: u64) -> bool {
    scan::window_owns_record_start(pos, win_start, win_end)
}

fn is_mtdna_chrom_norm(chrom: &str) -> bool {
    scan::is_mtdna_chrom_norm(chrom)
}

fn is_rdna_chrom_norm(chrom: &str, rdna_contigs: &HashSet<String>) -> bool {
    scan::is_rdna_chrom_norm(chrom, rdna_contigs)
}

fn scan_inline_qc_sample(
    bam_path: &str,
    header: &noodles::sam::Header,
    dense_maps: &HashMap<String, Arc<DenseMap>>,
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
        let Some(dense) = dense_maps.get(chrom_norm) else {
            continue;
        };

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
