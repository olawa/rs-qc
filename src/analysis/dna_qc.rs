use crate::analysis::bam_scan::scan_bam_stream;
use crate::analysis::report::write_summary_json;
use crate::io::bam::{for_each_aligned_block, normalized_reference_name};
use crate::models::normalize_chrom;
use anyhow::{bail, Context, Result};
use noodles::bam;
use noodles::sam;
use serde::Serialize;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

const BREADTH_THRESHOLDS: [u32; 5] = [1, 5, 10, 20, 30];

#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct DnaQcConfig {
    pub inputs: Vec<String>,
    pub output_prefix: String,
    pub mapq_threshold: u8,
    pub threads: usize,
    pub window_size: u32,
    pub targets_path: Option<String>,
    pub thresholds: Vec<u32>,
    pub callable_depth: u32,
    pub include_duplicates: bool,
    pub reference_fasta: Option<String>,
    pub annotation_path: Option<String>,
    pub low_cov_threshold: u32,
    pub snap_lowcov: u32,
    pub show_progress: bool,
}

#[derive(Debug, Clone, Default, Serialize, PartialEq)]
pub struct GcBiasBin {
    pub gc_fraction: f64,
    pub window_count: u64,
    pub mean_depth: f64,
    pub normalized_coverage: f64,
    pub standard_error: f64,
}

#[derive(Debug, Clone, Default, Serialize, PartialEq)]
pub struct GcBiasMetrics {
    pub bins: Vec<GcBiasBin>,
}

#[derive(Debug, Clone, Default)]
pub struct DnaQcMetrics {
    pub total_records: u64,
    pub mapped_records: u64,
    pub filtered_mapq_records: u64,
    pub total_aligned_bases: u64,
    pub total_reference_bases: u64,
    pub depth_hist: BTreeMap<u32, u64>,
    pub contigs: BTreeMap<String, DnaContigAccumulator>,
    pub windows: BTreeMap<String, DnaWindowAccumulator>,
    pub targets: BTreeMap<String, DnaTargetAccumulator>,
    pub target_bed: Option<String>,
    pub window_size: u32,
    pub gc_bias: Option<GcBiasMetrics>,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct DnaContigAccumulator {
    pub chrom: String,
    pub length: u64,
    pub aligned_bases: u64,
    pub depth_bases: u64,
    pub bases_ge_1x: u64,
    pub bases_ge_5x: u64,
    pub bases_ge_10x: u64,
    pub bases_ge_20x: u64,
    pub bases_ge_30x: u64,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct DnaWindowAccumulator {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub aligned_bases: u64,
    pub depth_bases: u64,
    pub bases_ge_1x: u64,
    pub bases_ge_5x: u64,
    pub bases_ge_10x: u64,
    pub bases_ge_20x: u64,
    pub bases_ge_30x: u64,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct DnaTargetAccumulator {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub aligned_bases: u64,
    pub depth_bases: u64,
    pub bases_ge_1x: u64,
    pub bases_ge_5x: u64,
    pub bases_ge_10x: u64,
    pub bases_ge_20x: u64,
    pub bases_ge_30x: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct DnaQcSummary {
    pub total_records: u64,
    pub mapped_records: u64,
    pub filtered_mapq_records: u64,
    pub total_aligned_bases: u64,
    pub total_reference_bases: u64,
    pub mean_depth: f64,
    pub median_depth: f64,
    pub breadth_1x: f64,
    pub breadth_5x: f64,
    pub breadth_10x: f64,
    pub breadth_20x: f64,
    pub breadth_30x: f64,
    pub depth_hist: BTreeMap<u32, u64>,
    pub contigs: Vec<DnaContigAccumulator>,
    pub windows: Vec<DnaWindowAccumulator>,
    pub targets: Vec<DnaTargetAccumulator>,
    pub target_bed: Option<String>,
    pub window_size: u32,
    pub gc_bias: Option<GcBiasMetrics>,
}

#[derive(Debug, Clone)]
struct TargetSpec {
    chrom: String,
    start: u64,
    end: u64,
    name: String,
}

#[derive(Debug, Clone)]
struct ContigState {
    chrom: String,
    length: u64,
    cursor: u64,
    depth: u32,
    end_heap: BinaryHeap<Reverse<u64>>,
    window_cursor: usize,
    target_cursor: usize,
}

impl DnaQcMetrics {
    pub fn summary(&self) -> DnaQcSummary {
        let mean_depth = if self.total_reference_bases == 0 {
            0.0
        } else {
            self.total_aligned_bases as f64 / self.total_reference_bases as f64
        };
        let median_depth = depth_quantile(&self.depth_hist, self.total_reference_bases, 0.5);
        let breadth = |threshold: u32| -> f64 {
            let bases = breadth_bases(&self.depth_hist, threshold);
            if self.total_reference_bases == 0 {
                0.0
            } else {
                bases as f64 * 100.0 / self.total_reference_bases as f64
            }
        };

        DnaQcSummary {
            total_records: self.total_records,
            mapped_records: self.mapped_records,
            filtered_mapq_records: self.filtered_mapq_records,
            total_aligned_bases: self.total_aligned_bases,
            total_reference_bases: self.total_reference_bases,
            mean_depth,
            median_depth,
            breadth_1x: breadth(BREADTH_THRESHOLDS[0]),
            breadth_5x: breadth(BREADTH_THRESHOLDS[1]),
            breadth_10x: breadth(BREADTH_THRESHOLDS[2]),
            breadth_20x: breadth(BREADTH_THRESHOLDS[3]),
            breadth_30x: breadth(BREADTH_THRESHOLDS[4]),
            depth_hist: self.depth_hist.clone(),
            contigs: self.contigs.values().cloned().collect(),
            windows: self.windows.values().cloned().collect(),
            targets: self.targets.values().cloned().collect(),
            target_bed: self.target_bed.clone(),
            window_size: self.window_size,
            gc_bias: self.gc_bias.clone(),
        }
    }

    fn merge(&mut self, other: DnaQcMetrics) {
        self.total_records += other.total_records;
        self.mapped_records += other.mapped_records;
        self.filtered_mapq_records += other.filtered_mapq_records;
        self.total_aligned_bases += other.total_aligned_bases;
        self.total_reference_bases += other.total_reference_bases;

        for (depth, count) in other.depth_hist {
            *self.depth_hist.entry(depth).or_insert(0) += count;
        }

        merge_accumulator_maps(&mut self.contigs, other.contigs);
        merge_accumulator_maps(&mut self.windows, other.windows);
        merge_accumulator_maps(&mut self.targets, other.targets);

        if self.gc_bias.is_none() {
            self.gc_bias = other.gc_bias;
        }
    }
}

impl DnaContigAccumulator {
    fn add_segment(&mut self, start: u64, end: u64, depth: u32) {
        if end <= start {
            return;
        }
        let len = end - start;
        self.depth_bases += len * u64::from(depth);
        if depth >= 1 {
            self.bases_ge_1x += len;
        }
        if depth >= 5 {
            self.bases_ge_5x += len;
        }
        if depth >= 10 {
            self.bases_ge_10x += len;
        }
        if depth >= 20 {
            self.bases_ge_20x += len;
        }
        if depth >= 30 {
            self.bases_ge_30x += len;
        }
    }
}

impl DnaWindowAccumulator {
    fn add_segment(&mut self, start: u64, end: u64, depth: u32) {
        if end <= start {
            return;
        }
        let len = end - start;
        self.depth_bases += len * u64::from(depth);
        self.aligned_bases += len * u64::from(depth > 0);
        if depth >= 1 {
            self.bases_ge_1x += len;
        }
        if depth >= 5 {
            self.bases_ge_5x += len;
        }
        if depth >= 10 {
            self.bases_ge_10x += len;
        }
        if depth >= 20 {
            self.bases_ge_20x += len;
        }
        if depth >= 30 {
            self.bases_ge_30x += len;
        }
    }
}

impl DnaTargetAccumulator {
    fn add_segment(&mut self, start: u64, end: u64, depth: u32) {
        if end <= start {
            return;
        }
        let len = end - start;
        self.depth_bases += len * u64::from(depth);
        self.aligned_bases += len * u64::from(depth > 0);
        if depth >= 1 {
            self.bases_ge_1x += len;
        }
        if depth >= 5 {
            self.bases_ge_5x += len;
        }
        if depth >= 10 {
            self.bases_ge_10x += len;
        }
        if depth >= 20 {
            self.bases_ge_20x += len;
        }
        if depth >= 30 {
            self.bases_ge_30x += len;
        }
    }
}

pub fn run_dna_qc(config: &DnaQcConfig) -> Result<DnaQcMetrics> {
    if config.inputs.is_empty() {
        bail!("at least one BAM input is required");
    }
    if config.window_size == 0 {
        bail!("window size must be greater than zero");
    }

    let mut metrics = DnaQcMetrics {
        target_bed: config.targets_path.clone(),
        window_size: config.window_size,
        ..DnaQcMetrics::default()
    };

    for input in &config.inputs {
        let file_metrics = scan_dna_file(input, config)
            .with_context(|| format!("failed while scanning DNA input {input}"))?;
        metrics.merge(file_metrics);
    }

    if let Some(ref ref_fasta) = config.reference_fasta {
        println!("Calculating GC Bias metrics using reference FASTA: {}", ref_fasta);
        let gc_bias = compute_gc_bias(ref_fasta, &metrics.windows)?;
        metrics.gc_bias = Some(gc_bias);
    }

    write_outputs(config, &metrics)?;
    Ok(metrics)
}

fn scan_dna_file(path: &str, config: &DnaQcConfig) -> Result<DnaQcMetrics> {
    if path.ends_with(".cram") {
        bail!("CRAM input is planned, but this build currently supports BAM for rs-qc dna");
    }

    let header = read_header(path)?;
    let state = DnaCoverageState::new(&header, config.window_size, config.targets_path.as_deref())?;
    let metrics = scan_bam_stream(
        path,
        &crate::analysis::bam_scan::BamScanConfig {
            threads: config.threads,
            show_progress: config.show_progress,
        },
        state,
        |state, header, record| {
            state.observe_record(header, record, 0);
        },
    )?;

    Ok(metrics.finish())
}

fn read_header(path: &str) -> Result<sam::Header> {
    let file = File::open(path).with_context(|| format!("could not open {path}"))?;
    let mut reader = bam::io::Reader::new(file);
    Ok(reader.read_header()?)
}

struct DnaCoverageState {
    metrics: DnaQcMetrics,
    reference_lengths: HashMap<String, u64>,
    windows_by_chrom: BTreeMap<String, Vec<String>>,
    targets_by_chrom: BTreeMap<String, Vec<String>>,
    current: Option<ContigState>,
}

impl DnaCoverageState {
    fn new(header: &sam::Header, window_size: u32, target_bed: Option<&str>) -> Result<Self> {
        let mut metrics = DnaQcMetrics {
            window_size,
            ..DnaQcMetrics::default()
        };
        let mut reference_lengths = HashMap::new();
        let mut windows_by_chrom: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut total_reference_bases = 0_u64;

        for (name, seq) in header.reference_sequences() {
            let chrom =
                normalize_chrom(String::from_utf8_lossy(name.as_ref()).as_ref()).into_owned();
            let length = seq.length().get() as u64;
            total_reference_bases += length;
            reference_lengths.insert(chrom.clone(), length);
            metrics.contigs.insert(
                chrom.clone(),
                DnaContigAccumulator {
                    chrom: chrom.clone(),
                    length,
                    ..DnaContigAccumulator::default()
                },
            );

            let mut windows = Vec::new();
            let mut start = 0_u64;
            while start < length {
                let end = (start + window_size as u64).min(length);
                let key = format!("{chrom}:{start}-{end}");
                metrics.windows.insert(
                    key.clone(),
                    DnaWindowAccumulator {
                        chrom: chrom.clone(),
                        start,
                        end,
                        ..DnaWindowAccumulator::default()
                    },
                );
                windows.push(key);
                start = end;
            }
            windows_by_chrom.insert(chrom, windows);
        }

        metrics.total_reference_bases = total_reference_bases;
        metrics.depth_hist.insert(0, total_reference_bases);

        let mut targets_by_chrom: BTreeMap<String, Vec<String>> = BTreeMap::new();
        if let Some(path) = target_bed {
            let targets = load_targets(path)?;
            metrics.target_bed = Some(path.to_string());
            for target in targets {
                let key = format!(
                    "{}:{}-{}:{}",
                    target.chrom, target.start, target.end, target.name
                );
                metrics.targets.insert(
                    key.clone(),
                    DnaTargetAccumulator {
                        chrom: target.chrom.clone(),
                        start: target.start,
                        end: target.end,
                        name: target.name,
                        ..DnaTargetAccumulator::default()
                    },
                );
                targets_by_chrom.entry(target.chrom).or_default().push(key);
            }
        }

        Ok(Self {
            metrics,
            reference_lengths,
            windows_by_chrom,
            targets_by_chrom,
            current: None,
        })
    }

    fn finish(mut self) -> DnaQcMetrics {
        self.flush_current_contig();
        self.metrics
    }

    fn observe_record(&mut self, header: &sam::Header, record: &bam::Record, mapq_threshold: u8) {
        self.metrics.total_records += 1;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            return;
        }
        self.metrics.mapped_records += 1;

        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < mapq_threshold {
            self.metrics.filtered_mapq_records += 1;
            return;
        }

        let Some(chrom) = normalized_reference_name(header, record) else {
            return;
        };
        let Some(contig_len) = self.reference_lengths.get(&chrom).copied() else {
            return;
        };
        self.ensure_current(&chrom, contig_len);

        for_each_aligned_block(record, |start, end| {
            self.metrics.total_aligned_bases += end - start;
            self.advance_to(start);
            self.push_interval(end);
        });
    }

    fn ensure_current(&mut self, chrom: &str, length: u64) {
        let needs_switch = self
            .current
            .as_ref()
            .map(|state| state.chrom != chrom)
            .unwrap_or(true);
        if needs_switch {
            self.flush_current_contig();
            self.current = Some(ContigState {
                chrom: chrom.to_string(),
                length,
                cursor: 0,
                depth: 0,
                end_heap: BinaryHeap::new(),
                window_cursor: 0,
                target_cursor: 0,
            });
        }
    }

    fn flush_current_contig(&mut self) {
        let Some(mut state) = self.current.take() else {
            return;
        };
        let length = state.length;
        self.advance_state(&mut state, length);
    }

    fn push_interval(&mut self, end: u64) {
        if let Some(state) = self.current.as_mut() {
            state.end_heap.push(Reverse(end));
            state.depth = state.depth.saturating_add(1);
        }
    }

    fn advance_to(&mut self, pos: u64) {
        let Some(mut state) = self.current.take() else {
            return;
        };
        self.advance_state(&mut state, pos);
        self.current = Some(state);
    }

    fn advance_state(&mut self, state: &mut ContigState, pos: u64) {
        while let Some(&Reverse(end)) = state.end_heap.peek() {
            if end > pos {
                break;
            }
            self.finalize_segment(
                &state.chrom,
                state.cursor,
                end,
                state.depth,
                &mut state.window_cursor,
                &mut state.target_cursor,
            );
            state.cursor = end;
            state.end_heap.pop();
            state.depth = state.depth.saturating_sub(1);
        }

        if pos > state.cursor {
            self.finalize_segment(
                &state.chrom,
                state.cursor,
                pos,
                state.depth,
                &mut state.window_cursor,
                &mut state.target_cursor,
            );
            state.cursor = pos;
        }
    }

    fn finalize_segment(
        &mut self,
        chrom: &str,
        start: u64,
        end: u64,
        depth: u32,
        window_cursor: &mut usize,
        target_cursor: &mut usize,
    ) {
        if end <= start {
            return;
        }

        let len = end - start;
        if depth == 0 {
            return;
        }
        if let Some(zero) = self.metrics.depth_hist.get_mut(&0) {
            *zero = zero.saturating_sub(len);
        }
        *self.metrics.depth_hist.entry(depth).or_insert(0) += len;

        if let Some(contig) = self.metrics.contigs.get_mut(chrom) {
            contig.add_segment(start, end, depth);
        }

        if let Some(keys) = self.windows_by_chrom.get(chrom) {
            let mut idx = *window_cursor;
            while idx < keys.len() {
                let key = &keys[idx];
                let Some(window) = self.metrics.windows.get_mut(key) else {
                    idx += 1;
                    continue;
                };
                if window.end <= start {
                    idx += 1;
                    continue;
                }
                if window.start >= end {
                    break;
                }
                let seg_start = start.max(window.start);
                let seg_end = end.min(window.end);
                if seg_end > seg_start {
                    window.add_segment(seg_start, seg_end, depth);
                }
                if window.end <= end {
                    idx += 1;
                } else {
                    break;
                }
            }
            *window_cursor = idx;
        }

        if let Some(keys) = self.targets_by_chrom.get(chrom) {
            let mut idx = *target_cursor;
            while idx < keys.len() {
                let key = &keys[idx];
                let Some(target) = self.metrics.targets.get_mut(key) else {
                    idx += 1;
                    continue;
                };
                if target.end <= start {
                    idx += 1;
                    continue;
                }
                if target.start >= end {
                    break;
                }
                let seg_start = start.max(target.start);
                let seg_end = end.min(target.end);
                if seg_end > seg_start {
                    target.add_segment(seg_start, seg_end, depth);
                }
                if target.end <= end {
                    idx += 1;
                } else {
                    break;
                }
            }
            *target_cursor = idx;
        }
    }
}

fn load_targets(path: &str) -> Result<Vec<TargetSpec>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open DNA targets BED file: {}", path))?;
    let reader = BufReader::new(file);
    let mut targets = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            bail!("Invalid BED line {} in {}", line_no + 1, path);
        }
        let chrom = normalize_chrom(fields[0]);
        let start: u64 = fields[1]
            .parse()
            .with_context(|| format!("Invalid BED start at line {}", line_no + 1))?;
        let end: u64 = fields[2]
            .parse()
            .with_context(|| format!("Invalid BED end at line {}", line_no + 1))?;
        let name = fields
            .get(3)
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, start, end));
        if start < end {
            targets.push(TargetSpec {
                chrom: chrom.into_owned(),
                start,
                end,
                name,
            });
        }
    }

    targets.sort_by(|a, b| a.chrom.cmp(&b.chrom).then_with(|| a.start.cmp(&b.start)));
    Ok(targets)
}

fn merge_accumulator_maps<T: Clone>(dst: &mut BTreeMap<String, T>, src: BTreeMap<String, T>)
where
    T: MergeableAccumulator,
{
    for (key, value) in src {
        dst.entry(key)
            .and_modify(|existing| existing.merge_from(&value))
            .or_insert(value);
    }
}

trait MergeableAccumulator {
    fn merge_from(&mut self, other: &Self);
}

impl MergeableAccumulator for DnaContigAccumulator {
    fn merge_from(&mut self, other: &Self) {
        self.length = self.length.max(other.length);
        self.aligned_bases += other.aligned_bases;
        self.depth_bases += other.depth_bases;
        self.bases_ge_1x += other.bases_ge_1x;
        self.bases_ge_5x += other.bases_ge_5x;
        self.bases_ge_10x += other.bases_ge_10x;
        self.bases_ge_20x += other.bases_ge_20x;
        self.bases_ge_30x += other.bases_ge_30x;
    }
}

impl MergeableAccumulator for DnaWindowAccumulator {
    fn merge_from(&mut self, other: &Self) {
        self.aligned_bases += other.aligned_bases;
        self.depth_bases += other.depth_bases;
        self.bases_ge_1x += other.bases_ge_1x;
        self.bases_ge_5x += other.bases_ge_5x;
        self.bases_ge_10x += other.bases_ge_10x;
        self.bases_ge_20x += other.bases_ge_20x;
        self.bases_ge_30x += other.bases_ge_30x;
    }
}

impl MergeableAccumulator for DnaTargetAccumulator {
    fn merge_from(&mut self, other: &Self) {
        self.aligned_bases += other.aligned_bases;
        self.depth_bases += other.depth_bases;
        self.bases_ge_1x += other.bases_ge_1x;
        self.bases_ge_5x += other.bases_ge_5x;
        self.bases_ge_10x += other.bases_ge_10x;
        self.bases_ge_20x += other.bases_ge_20x;
        self.bases_ge_30x += other.bases_ge_30x;
    }
}

fn breadth_bases(depth_hist: &BTreeMap<u32, u64>, threshold: u32) -> u64 {
    depth_hist
        .iter()
        .filter(|(depth, _)| **depth >= threshold)
        .map(|(_, bases)| *bases)
        .sum()
}

fn depth_quantile(depth_hist: &BTreeMap<u32, u64>, total_bases: u64, q: f64) -> f64 {
    if total_bases == 0 {
        return 0.0;
    }
    let target = ((total_bases as f64 - 1.0) * q.clamp(0.0, 1.0)).round() as u64 + 1;
    let mut seen = 0_u64;
    for (&depth, &bases) in depth_hist {
        seen += bases;
        if seen >= target {
            return depth as f64;
        }
    }
    0.0
}

fn write_outputs(config: &DnaQcConfig, metrics: &DnaQcMetrics) -> Result<()> {
    let summary = metrics.summary();
    write_summary(config, &summary)?;
    write_depth_hist(config, &summary)?;
    write_contigs(config, &summary)?;
    write_windows(config, &summary)?;
    write_targets(config, &summary)?;
    if let Some(ref gc) = metrics.gc_bias {
        write_gc_bias(config, gc)?;
    }
    write_summary_json(
        &format!("{}.dna.summary.json", config.output_prefix),
        "dna",
        &config.output_prefix,
        &summary,
    )?;
    Ok(())
}

fn write_summary(config: &DnaQcConfig, summary: &DnaQcSummary) -> Result<()> {
    let mut out = File::create(format!("{}.dna.summary.tsv", config.output_prefix))?;
    writeln!(out, "metric\tvalue")?;
    writeln!(out, "total_records\t{}", summary.total_records)?;
    writeln!(out, "mapped_records\t{}", summary.mapped_records)?;
    writeln!(
        out,
        "filtered_mapq_records\t{}",
        summary.filtered_mapq_records
    )?;
    writeln!(out, "total_aligned_bases\t{}", summary.total_aligned_bases)?;
    writeln!(
        out,
        "total_reference_bases\t{}",
        summary.total_reference_bases
    )?;
    writeln!(out, "mean_depth\t{:.6}", summary.mean_depth)?;
    writeln!(out, "median_depth\t{:.6}", summary.median_depth)?;
    writeln!(out, "breadth_1x\t{:.4}", summary.breadth_1x)?;
    writeln!(out, "breadth_5x\t{:.4}", summary.breadth_5x)?;
    writeln!(out, "breadth_10x\t{:.4}", summary.breadth_10x)?;
    writeln!(out, "breadth_20x\t{:.4}", summary.breadth_20x)?;
    writeln!(out, "breadth_30x\t{:.4}", summary.breadth_30x)?;
    writeln!(out, "window_size\t{}", summary.window_size)?;
    if let Some(target_bed) = &summary.target_bed {
        writeln!(out, "target_bed\t{}", target_bed)?;
    }
    Ok(())
}

fn write_depth_hist(config: &DnaQcConfig, summary: &DnaQcSummary) -> Result<()> {
    let mut out = File::create(format!("{}.dna.depth_hist.tsv", config.output_prefix))?;
    writeln!(out, "depth\tbases\tfraction")?;
    let total = summary.total_reference_bases.max(1);
    for (depth, bases) in &summary.depth_hist {
        writeln!(
            out,
            "{}\t{}\t{:.8}",
            depth,
            bases,
            *bases as f64 / total as f64
        )?;
    }
    Ok(())
}

fn write_contigs(config: &DnaQcConfig, summary: &DnaQcSummary) -> Result<()> {
    let mut out = File::create(format!("{}.dna.contigs.tsv", config.output_prefix))?;
    writeln!(
        out,
        "chrom\tlength\taligned_bases\tdepth_bases\tmean_depth\tbreadth_1x\tbreadth_5x\tbreadth_10x\tbreadth_20x\tbreadth_30x"
    )?;
    for contig in &summary.contigs {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.6}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            contig.chrom,
            contig.length,
            contig.aligned_bases,
            contig.depth_bases,
            if contig.length == 0 {
                0.0
            } else {
                contig.depth_bases as f64 / contig.length as f64
            },
            pct(contig.bases_ge_1x, contig.length),
            pct(contig.bases_ge_5x, contig.length),
            pct(contig.bases_ge_10x, contig.length),
            pct(contig.bases_ge_20x, contig.length),
            pct(contig.bases_ge_30x, contig.length),
        )?;
    }
    Ok(())
}

fn write_windows(config: &DnaQcConfig, summary: &DnaQcSummary) -> Result<()> {
    let mut out = File::create(format!("{}.dna.windows.tsv", config.output_prefix))?;
    writeln!(
        out,
        "chrom\tstart\tend\taligned_bases\tdepth_bases\tmean_depth\tbreadth_1x\tbreadth_5x\tbreadth_10x\tbreadth_20x\tbreadth_30x"
    )?;
    for window in &summary.windows {
        let len = window.end.saturating_sub(window.start).max(1);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            window.chrom,
            window.start,
            window.end,
            window.aligned_bases,
            window.depth_bases,
            window.depth_bases as f64 / len as f64,
            pct(window.bases_ge_1x, len),
            pct(window.bases_ge_5x, len),
            pct(window.bases_ge_10x, len),
            pct(window.bases_ge_20x, len),
            pct(window.bases_ge_30x, len),
        )?;
    }
    Ok(())
}

fn write_targets(config: &DnaQcConfig, summary: &DnaQcSummary) -> Result<()> {
    let mut out = File::create(format!("{}.dna.targets.tsv", config.output_prefix))?;
    writeln!(
        out,
        "chrom\tstart\tend\tname\taligned_bases\tdepth_bases\tmean_depth\tbreadth_1x\tbreadth_5x\tbreadth_10x\tbreadth_20x\tbreadth_30x"
    )?;
    for target in &summary.targets {
        let len = target.end.saturating_sub(target.start).max(1);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            target.chrom,
            target.start,
            target.end,
            target.name,
            target.aligned_bases,
            target.depth_bases,
            target.depth_bases as f64 / len as f64,
            pct(target.bases_ge_1x, len),
            pct(target.bases_ge_5x, len),
            pct(target.bases_ge_10x, len),
            pct(target.bases_ge_20x, len),
            pct(target.bases_ge_30x, len),
        )?;
    }
    Ok(())
}

fn pct(n: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        n as f64 * 100.0 / total as f64
    }
}

pub fn compute_gc_bias(
    reference_fasta: &str,
    windows: &BTreeMap<String, DnaWindowAccumulator>,
) -> Result<GcBiasMetrics> {
    let mut reader = noodles::fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_fasta)
        .with_context(|| format!("could not open indexed FASTA {} (make sure the .fai index file exists)", reference_fasta))?;

    let mut gc_bin_depths: Vec<Vec<f64>> = vec![Vec::new(); 101];

    for window in windows.values() {
        let len = window.end - window.start;
        if len == 0 {
            continue;
        }

        if let Some(seq) = get_fasta_sequence(&mut reader, &window.chrom, window.start, window.end) {
            let mut gc_count = 0;
            let mut at_count = 0;
            for &b in &seq {
                match b.to_ascii_uppercase() {
                    b'G' | b'C' => gc_count += 1,
                    b'A' | b'T' | b'U' => at_count += 1,
                    _ => {}
                }
            }
            let total = gc_count + at_count;
            if total > 0 {
                let gc_frac = gc_count as f64 / total as f64;
                let bin_idx = (gc_frac * 100.0).round() as usize;
                let bin_idx = bin_idx.min(100);

                let mean_depth = window.depth_bases as f64 / len as f64;
                gc_bin_depths[bin_idx].push(mean_depth);
            }
        }
    }

    let all_depths: Vec<f64> = gc_bin_depths.iter().flatten().copied().collect();
    let total_windows = all_depths.len() as u64;
    let global_mean_depth = if total_windows == 0 {
        0.0
    } else {
        all_depths.iter().sum::<f64>() / total_windows as f64
    };

    let mut bins = Vec::with_capacity(101);
    for (gc_idx, depths) in gc_bin_depths.iter().enumerate() {
        let count = depths.len() as u64;
        let (mean_depth, normalized_coverage, standard_error) = if count == 0 {
            (0.0, 0.0, 0.0)
        } else {
            let sum: f64 = depths.iter().sum();
            let mean = sum / count as f64;

            let norm_cov = if global_mean_depth > 0.0 {
                mean / global_mean_depth
            } else {
                0.0
            };

            let se = if count > 1 {
                let variance_sum: f64 = depths.iter().map(|&d| (d - mean).powi(2)).sum();
                let variance = variance_sum / (count - 1) as f64;
                variance.sqrt() / (count as f64).sqrt()
            } else {
                0.0
            };

            (mean, norm_cov, se)
        };

        bins.push(GcBiasBin {
            gc_fraction: gc_idx as f64 / 100.0,
            window_count: count,
            mean_depth,
            normalized_coverage,
            standard_error,
        });
    }

    // Run premium validation check outputs in the console
    let warnings = run_gc_bias_validation(&bins);
    if !warnings.is_empty() {
        println!("  [WARNING] DNA GC Bias Validation warnings found:");
        for w in &warnings {
            println!("    - {}", w);
        }
    } else {
        println!("  - DNA GC Bias Validation: PASS (no extreme GC or AT coverage drop-offs detected)");
    }

    Ok(GcBiasMetrics { bins })
}

fn get_fasta_sequence<R>(
    reader: &mut noodles::fasta::io::indexed_reader::IndexedReader<R>,
    chrom: &str,
    start: u64,
    end: u64,
) -> Option<Vec<u8>>
where
    R: std::io::BufRead + std::io::Seek,
{
    let trimmed = chrom.trim();
    let lower = trimmed.to_ascii_lowercase();
    let no_chr = lower.strip_prefix("chr").unwrap_or(&lower);
    let canonical_with_chr = format!("chr{}", no_chr);
    let canonical_without_chr = no_chr.to_string();

    let aliases = vec![
        trimmed.to_string(),
        canonical_with_chr,
        canonical_without_chr,
    ];

    for alias in aliases {
        let region_str = format!("{}:{}-{}", alias, start + 1, end);
        if let Ok(region) = region_str.parse() {
            if let Ok(rec) = reader.query(&region) {
                return Some(rec.sequence().as_ref().to_vec());
            }
        }
    }
    None
}

pub fn run_gc_bias_validation(bins: &[GcBiasBin]) -> Vec<String> {
    let mut warnings = Vec::new();

    // Check AT-rich regions (e.g. GC bin 20%-25%)
    let mut at_sum = 0.0;
    let mut at_count = 0;
    for bin in bins {
        let pct = (bin.gc_fraction * 100.0).round() as u32;
        if (20..=25).contains(&pct) && bin.window_count > 0 {
            at_sum += bin.normalized_coverage;
            at_count += 1;
        }
    }
    if at_count > 0 {
        let at_avg = at_sum / at_count as f64;
        if at_avg < 0.5 {
            warnings.push(format!("AT-rich regions (20-25% GC) have extremely depleted coverage: {:.2}x normalized", at_avg));
        } else if at_avg > 1.8 {
            warnings.push(format!("AT-rich regions (20-25% GC) have extreme coverage spikes: {:.2}x normalized", at_avg));
        }
    }

    // Check GC-rich regions (e.g. GC bin 60%-65%)
    let mut gc_sum = 0.0;
    let mut gc_count = 0;
    for bin in bins {
        let pct = (bin.gc_fraction * 100.0).round() as u32;
        if (60..=65).contains(&pct) && bin.window_count > 0 {
            gc_sum += bin.normalized_coverage;
            gc_count += 1;
        }
    }
    if gc_count > 0 {
        let gc_avg = gc_sum / gc_count as f64;
        if gc_avg < 0.5 {
            warnings.push(format!("GC-rich regions (60-65% GC) have extremely depleted coverage: {:.2}x normalized", gc_avg));
        } else if gc_avg > 1.8 {
            warnings.push(format!("GC-rich regions (60-65% GC) have extreme coverage spikes: {:.2}x normalized", gc_avg));
        }
    }

    // Check overall standard error mean across bins that contain windows
    let mut total_se = 0.0;
    let mut bins_with_windows = 0;
    for bin in bins {
        if bin.window_count > 1 {
            total_se += bin.standard_error;
            bins_with_windows += 1;
        }
    }
    if bins_with_windows > 0 {
        let avg_se = total_se / bins_with_windows as f64;
        if avg_se > 0.15 {
            warnings.push(format!("High average standard error across GC bins ({:.4}): depth estimates within GC groups are highly variable/noisy", avg_se));
        }
    }

    warnings
}

fn write_gc_bias(config: &DnaQcConfig, metrics: &GcBiasMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.dna.gc_bias.tsv", config.output_prefix))?;
    writeln!(
        out,
        "gc_fraction\twindow_count\tmean_depth\tnormalized_coverage\tstandard_error"
    )?;
    for bin in &metrics.bins {
        writeln!(
            out,
            "{:.2}\t{}\t{:.6}\t{:.6}\t{:.6}",
            bin.gc_fraction,
            bin.window_count,
            bin.mean_depth,
            bin.normalized_coverage,
            bin.standard_error
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_bias_validation() {
        // 1. Pass case
        let mut normal_bins = Vec::new();
        for i in 0..=100 {
            normal_bins.push(GcBiasBin {
                gc_fraction: i as f64 / 100.0,
                window_count: 5,
                mean_depth: 30.0,
                normalized_coverage: 1.0,
                standard_error: 0.02,
            });
        }
        let pass_warnings = run_gc_bias_validation(&normal_bins);
        assert!(pass_warnings.is_empty(), "Expected no warnings for uniform coverage");

        // 2. AT depletion case (GC bin 20%-25% having low normalized coverage)
        let mut at_depleted_bins = normal_bins.clone();
        for bin in &mut at_depleted_bins {
            let pct = (bin.gc_fraction * 100.0).round() as u32;
            if (20..=25).contains(&pct) {
                bin.normalized_coverage = 0.4;
            }
        }
        let at_warnings = run_gc_bias_validation(&at_depleted_bins);
        assert!(!at_warnings.is_empty());
        assert!(at_warnings[0].contains("AT-rich"));

        // 3. GC depletion case (GC bin 60%-65% having low normalized coverage)
        let mut gc_depleted_bins = normal_bins.clone();
        for bin in &mut gc_depleted_bins {
            let pct = (bin.gc_fraction * 100.0).round() as u32;
            if (60..=65).contains(&pct) {
                bin.normalized_coverage = 0.3;
            }
        }
        let gc_warnings = run_gc_bias_validation(&gc_depleted_bins);
        assert!(!gc_warnings.is_empty());
        assert!(gc_warnings[0].contains("GC-rich"));

        // 4. High standard error case
        let mut noisy_bins = normal_bins.clone();
        for bin in &mut noisy_bins {
            bin.standard_error = 0.25;
        }
        let noisy_warnings = run_gc_bias_validation(&noisy_bins);
        assert!(!noisy_warnings.is_empty());
        assert!(noisy_warnings[0].contains("standard error"));
    }
}
