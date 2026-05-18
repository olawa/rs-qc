use crate::analysis::bam_scan::{scan_bam_stream, BamScanConfig};
use crate::analysis::report::write_summary_json;
use crate::io::bam::ReferenceNames;
use anyhow::{bail, Context, Result};
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

const ACCURACY_SCALE: f64 = 10_000.0;

#[derive(Clone, Debug)]
pub struct AlignmentQcConfig {
    pub inputs: Vec<String>,
    pub output_prefix: String,
    pub mapq_threshold: u8,
    pub threads: usize,
    pub show_progress: bool,
}

#[derive(Clone, Debug, Default, Serialize, PartialEq)]
pub struct InsertSizeMetrics {
    pub total_inserts: u64,
    pub mono_nucleosomal_peak: Option<u32>,
    pub di_nucleosomal_peak: Option<u32>,
    pub mono_nucleosomal_count: u64,
    pub di_nucleosomal_count: u64,
    pub sub_nucleosomal_count: u64,
    pub short_cfdna_count: u64,
    pub mono_cfdna_count: u64,
    pub cfdna_ratio: Option<f64>,
    pub mono_di_ratio: Option<f64>,
    pub short_fraction: f64,
    pub mono_fraction: f64,
    pub di_fraction: f64,
    pub large_fragment_count: u64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct AlignmentQcMetrics {
    pub total_records: u64,
    pub mapped_records: u64,
    pub unmapped_records: u64,
    pub primary_records: u64,
    pub secondary_records: u64,
    pub supplementary_records: u64,
    pub duplicate_records: u64,
    pub qc_fail_records: u64,
    pub paired_records: u64,
    pub properly_paired_records: u64,
    pub singleton_records: u64,
    pub discordant_records: u64,
    pub orphan_records: u64,
    pub mapq_filtered_records: u64,
    pub total_query_bases: u64,
    pub total_aligned_match_bases: u64,
    pub soft_clipped_bases: u64,
    pub hard_clipped_bases: u64,
    pub mapq_hist: BTreeMap<u8, u64>,
    pub read_length_hist: BTreeMap<usize, u64>,
    pub insert_size_hist: BTreeMap<u32, u64>,
    pub cigar_op_bases: BTreeMap<String, u64>,
    pub per_contig: BTreeMap<String, ContigAlignmentMetrics>,
    pub accuracy: AccuracyMetrics,
    pub insert_size_metrics: Option<InsertSizeMetrics>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ContigAlignmentMetrics {
    pub records: u64,
    pub mapped_records: u64,
    pub primary_mapped_records: u64,
    pub bases: u64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct AccuracyMetrics {
    pub records_with_de: u64,
    pub accuracy_sum: f64,
    pub accuracy_hist: BTreeMap<u32, u64>,
}

impl AccuracyMetrics {
    pub fn mean(&self) -> Option<f64> {
        if self.records_with_de == 0 {
            None
        } else {
            Some(self.accuracy_sum / self.records_with_de as f64)
        }
    }

    pub fn median(&self) -> Option<f64> {
        self.quantile(0.5)
    }

    pub fn mode(&self) -> Option<f64> {
        self.accuracy_hist
            .iter()
            .max_by(|(left_bin, left_count), (right_bin, right_count)| {
                left_count
                    .cmp(right_count)
                    .then_with(|| right_bin.cmp(left_bin))
            })
            .map(|(&bin, _)| bin as f64 / ACCURACY_SCALE)
    }

    pub fn quantile(&self, q: f64) -> Option<f64> {
        if self.records_with_de == 0 {
            return None;
        }
        let target = ((self.records_with_de as f64 - 1.0) * q.clamp(0.0, 1.0)).round() as u64 + 1;
        let mut seen = 0_u64;
        for (&bin, &count) in &self.accuracy_hist {
            seen += count;
            if seen >= target {
                return Some(bin as f64 / ACCURACY_SCALE);
            }
        }
        None
    }

    fn observe(&mut self, de: f64) {
        let accuracy = (1.0 - de).clamp(0.0, 1.0);
        let bin = (accuracy * ACCURACY_SCALE).round() as u32;
        self.records_with_de += 1;
        self.accuracy_sum += accuracy;
        *self.accuracy_hist.entry(bin).or_insert(0) += 1;
    }
}

impl AlignmentQcMetrics {
    fn observe_record(
        &mut self,
        _header: &sam::Header,
        record: &bam::Record,
        mapq_threshold: u8,
        ref_names: &ReferenceNames,
    ) {
        self.total_records += 1;

        let flags = record.flags();
        if flags.is_unmapped() {
            self.unmapped_records += 1;
        } else {
            self.mapped_records += 1;
        }
        if flags.is_secondary() {
            self.secondary_records += 1;
        }
        if flags.is_supplementary() {
            self.supplementary_records += 1;
        }
        if !flags.is_secondary() && !flags.is_supplementary() {
            self.primary_records += 1;
        }
        if flags.is_duplicate() {
            self.duplicate_records += 1;
        }
        if flags.is_qc_fail() {
            self.qc_fail_records += 1;
        }
        if flags.is_segmented() {
            self.paired_records += 1;
            if flags.is_properly_segmented() {
                self.properly_paired_records += 1;
            }
            if flags.is_mate_unmapped() && !flags.is_unmapped() {
                self.singleton_records += 1;
            }
            if flags.is_unmapped() && !flags.is_mate_unmapped() {
                self.orphan_records += 1;
            }
            if !flags.is_properly_segmented() && !flags.is_unmapped() && !flags.is_mate_unmapped() {
                self.discordant_records += 1;
            }
        }

        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        *self.mapq_hist.entry(mapq).or_insert(0) += 1;
        if mapq < mapq_threshold {
            self.mapq_filtered_records += 1;
        }

        let read_len = record.sequence().len();
        self.total_query_bases += read_len as u64;
        *self.read_length_hist.entry(read_len).or_insert(0) += 1;

        let tlen = record.template_length();
        if flags.is_segmented() && tlen > 0 {
            *self.insert_size_hist.entry(tlen as u32).or_insert(0) += 1;
        }

        let mut aligned_match_bases = 0_u64;
        let mut reference_consuming_bases = 0_u64;
        for op_res in record.cigar().iter() {
            if let Ok(op) = op_res {
                let len = op.len() as u64;
                let kind = op.kind();
                let label = cigar_label(kind);

                // Use a more efficient entry pattern
                match self.cigar_op_bases.get_mut(label) {
                    Some(count) => *count += len,
                    None => {
                        self.cigar_op_bases.insert(label.to_string(), len);
                    }
                }

                match kind {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        aligned_match_bases += len;
                        reference_consuming_bases += len;
                    }
                    Kind::Deletion | Kind::Skip => {
                        reference_consuming_bases += len;
                    }
                    Kind::SoftClip => self.soft_clipped_bases += len,
                    Kind::HardClip => self.hard_clipped_bases += len,
                    _ => {}
                }
            }
        }
        self.total_aligned_match_bases += aligned_match_bases;

        let contig = ref_names
            .raw(record)
            .map(|s| s.to_string())
            .unwrap_or_else(|| {
                if flags.is_unmapped() {
                    "*".to_string()
                } else {
                    "unknown".to_string()
                }
            });
        let contig_metrics = self.per_contig.entry(contig).or_default();
        contig_metrics.records += 1;
        if !flags.is_unmapped() {
            contig_metrics.mapped_records += 1;
            contig_metrics.bases += reference_consuming_bases;
            if !flags.is_secondary() && !flags.is_supplementary() {
                contig_metrics.primary_mapped_records += 1;
            }
        }

        if let Some(de) = de_tag(record) {
            self.accuracy.observe(de as f64);
        }
    }
}

pub fn run_alignment_qc(config: &AlignmentQcConfig) -> Result<AlignmentQcMetrics> {
    if config.inputs.is_empty() {
        bail!("at least one BAM input is required");
    }

    let mut metrics = AlignmentQcMetrics::default();
    for input in &config.inputs {
        metrics = scan_alignment(input, config, metrics)
            .with_context(|| format!("failed while scanning alignment input {input}"))?;
    }

    metrics.insert_size_metrics = calculate_insert_size_metrics(&metrics.insert_size_hist);

    write_outputs(config, &metrics)?;
    Ok(metrics)
}

fn scan_alignment(
    path: &str,
    config: &AlignmentQcConfig,
    metrics: AlignmentQcMetrics,
) -> Result<AlignmentQcMetrics> {
    if path.ends_with(".cram") {
        bail!("CRAM input is planned, but this build currently supports BAM for rs-qc align");
    }

    let mut ref_names = None;
    scan_bam_stream(
        path,
        &BamScanConfig {
            threads: config.threads,
            show_progress: config.show_progress,
        },
        metrics,
        |state, header, record| {
            if ref_names.is_none() {
                ref_names = Some(ReferenceNames::new(header));
            }
            state.observe_record(
                header,
                record,
                config.mapq_threshold,
                ref_names.as_ref().unwrap(),
            );
        },
    )
}

fn write_outputs(config: &AlignmentQcConfig, metrics: &AlignmentQcMetrics) -> Result<()> {
    write_summary(config, metrics)?;
    write_histogram(
        &format!("{}.align.mapq.tsv", config.output_prefix),
        "mapq",
        "records",
        metrics.mapq_hist.iter().map(|(&k, &v)| (k as u64, v)),
    )?;
    write_histogram(
        &format!("{}.align.read_length.tsv", config.output_prefix),
        "read_length",
        "records",
        metrics
            .read_length_hist
            .iter()
            .map(|(&k, &v)| (k as u64, v)),
    )?;
    write_histogram(
        &format!("{}.align.insert_size.tsv", config.output_prefix),
        "insert_size",
        "records",
        metrics
            .insert_size_hist
            .iter()
            .map(|(&k, &v)| (k as u64, v)),
    )?;
    write_cigar(config, metrics)?;
    write_contigs(config, metrics)?;
    write_accuracy_hist(config, metrics)?;
    if let Some(ref m) = metrics.insert_size_metrics {
        write_insert_size_metrics(config, m)?;
    }
    write_summary_json(
        &format!("{}.align.summary.json", config.output_prefix),
        "align",
        &config.output_prefix,
        metrics,
    )?;
    Ok(())
}

fn write_summary(config: &AlignmentQcConfig, metrics: &AlignmentQcMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.align.summary.txt", config.output_prefix))?;
    writeln!(out, "total_records\t{}", metrics.total_records)?;
    writeln!(out, "mapped_records\t{}", metrics.mapped_records)?;
    writeln!(out, "unmapped_records\t{}", metrics.unmapped_records)?;
    writeln!(out, "primary_records\t{}", metrics.primary_records)?;
    writeln!(out, "secondary_records\t{}", metrics.secondary_records)?;
    writeln!(
        out,
        "supplementary_records\t{}",
        metrics.supplementary_records
    )?;
    writeln!(out, "duplicate_records\t{}", metrics.duplicate_records)?;
    writeln!(out, "qc_fail_records\t{}", metrics.qc_fail_records)?;
    writeln!(out, "paired_records\t{}", metrics.paired_records)?;
    writeln!(
        out,
        "properly_paired_records\t{}",
        metrics.properly_paired_records
    )?;
    writeln!(out, "singleton_records\t{}", metrics.singleton_records)?;
    writeln!(out, "discordant_records\t{}", metrics.discordant_records)?;
    writeln!(out, "orphan_records\t{}", metrics.orphan_records)?;
    writeln!(
        out,
        "mapq_filtered_records_lt_{}\t{}",
        config.mapq_threshold, metrics.mapq_filtered_records
    )?;
    writeln!(out, "total_query_bases\t{}", metrics.total_query_bases)?;
    writeln!(
        out,
        "total_aligned_match_bases\t{}",
        metrics.total_aligned_match_bases
    )?;
    writeln!(out, "soft_clipped_bases\t{}", metrics.soft_clipped_bases)?;
    writeln!(out, "hard_clipped_bases\t{}", metrics.hard_clipped_bases)?;
    writeln!(
        out,
        "de_accuracy_records\t{}",
        metrics.accuracy.records_with_de
    )?;
    writeln!(
        out,
        "de_accuracy_mean\t{}",
        fmt_optional(metrics.accuracy.mean())
    )?;
    writeln!(
        out,
        "de_accuracy_median\t{}",
        fmt_optional(metrics.accuracy.median())
    )?;
    writeln!(
        out,
        "de_accuracy_mode\t{}",
        fmt_optional(metrics.accuracy.mode())
    )?;
    writeln!(
        out,
        "de_accuracy_p05\t{}",
        fmt_optional(metrics.accuracy.quantile(0.05))
    )?;
    writeln!(
        out,
        "de_accuracy_p95\t{}",
        fmt_optional(metrics.accuracy.quantile(0.95))
    )?;
    Ok(())
}

fn write_cigar(config: &AlignmentQcConfig, metrics: &AlignmentQcMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.align.cigar.tsv", config.output_prefix))?;
    writeln!(out, "operation\tbases")?;
    for (op, bases) in &metrics.cigar_op_bases {
        writeln!(out, "{op}\t{bases}")?;
    }
    Ok(())
}

fn write_contigs(config: &AlignmentQcConfig, metrics: &AlignmentQcMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.align.contigs.tsv", config.output_prefix))?;
    writeln!(
        out,
        "contig\trecords\tmapped_records\tprimary_mapped_records\tbases"
    )?;
    for (contig, m) in &metrics.per_contig {
        writeln!(
            out,
            "{contig}\t{}\t{}\t{}\t{}",
            m.records, m.mapped_records, m.primary_mapped_records, m.bases
        )?;
    }
    Ok(())
}

fn write_accuracy_hist(config: &AlignmentQcConfig, metrics: &AlignmentQcMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.align.de_accuracy.tsv", config.output_prefix))?;
    writeln!(out, "accuracy_bin\tde_bin\trecords\tfraction")?;
    let total = metrics.accuracy.records_with_de.max(1);
    for (&bin, &count) in &metrics.accuracy.accuracy_hist {
        let accuracy = bin as f64 / ACCURACY_SCALE;
        let de = 1.0 - accuracy;
        writeln!(
            out,
            "{accuracy:.4}\t{de:.4}\t{count}\t{:.8}",
            count as f64 / total as f64
        )?;
    }
    Ok(())
}

fn write_histogram<I>(path: &str, key: &str, value: &str, rows: I) -> Result<()>
where
    I: IntoIterator<Item = (u64, u64)>,
{
    let mut out = File::create(path)?;
    writeln!(out, "{key}\t{value}")?;
    for (k, v) in rows {
        writeln!(out, "{k}\t{v}")?;
    }
    Ok(())
}

fn de_tag(record: &bam::Record) -> Option<f32> {
    let tag = Tag::new(b'd', b'e');
    match record.data().get(&tag)?.ok()? {
        Value::Float(v) => Some(v),
        Value::Int8(v) => Some(v as f32),
        Value::UInt8(v) => Some(v as f32),
        Value::Int16(v) => Some(v as f32),
        Value::UInt16(v) => Some(v as f32),
        Value::Int32(v) => Some(v as f32),
        Value::UInt32(v) => Some(v as f32),
        _ => None,
    }
}

fn cigar_label(kind: Kind) -> &'static str {
    match kind {
        Kind::Match => "M",
        Kind::Insertion => "I",
        Kind::Deletion => "D",
        Kind::Skip => "N",
        Kind::SoftClip => "S",
        Kind::HardClip => "H",
        Kind::Pad => "P",
        Kind::SequenceMatch => "=",
        Kind::SequenceMismatch => "X",
    }
}

fn fmt_optional(value: Option<f64>) -> String {
    value
        .map(|v| format!("{v:.6}"))
        .unwrap_or_else(|| "NA".to_string())
}

pub fn calculate_insert_size_metrics(hist: &BTreeMap<u32, u64>) -> Option<InsertSizeMetrics> {
    let total_inserts: u64 = hist.values().sum();
    if total_inserts == 0 {
        return None;
    }

    let mut sub_nucleosomal_count = 0_u64;
    let mut mono_nucleosomal_count = 0_u64;
    let mut di_nucleosomal_count = 0_u64;
    let mut short_cfdna_count = 0_u64;
    let mut mono_cfdna_count = 0_u64;
    let mut large_fragment_count = 0_u64;

    for (&size, &count) in hist {
        if size >= 30 && size < 143 {
            sub_nucleosomal_count += count;
        }
        if size >= 143 && size <= 220 {
            mono_nucleosomal_count += count;
        }
        if size >= 320 && size <= 480 {
            di_nucleosomal_count += count;
        }
        if size >= 100 && size <= 150 {
            short_cfdna_count += count;
        }
        if size >= 151 && size <= 220 {
            mono_cfdna_count += count;
        }
        if size > 1000 {
            large_fragment_count += count;
        }
    }

    // Helper for smoothing counts to identify peaks robustly (window size 7: -3..=3)
    let get_smoothed_count = |x: u32| -> f64 {
        let mut sum = 0.0;
        let mut count = 0;
        for offset in -3..=3 {
            let val = (x as i32 + offset) as u32;
            if let Some(&c) = hist.get(&val) {
                sum += c as f64;
            }
            count += 1;
        }
        sum / count as f64
    };

    // Find mono-nucleosomal peak (140-200 bp range)
    let mut mono_nucleosomal_peak = None;
    let mut max_mono_val = -1.0;
    for x in 140..=200 {
        let val = get_smoothed_count(x);
        if val > max_mono_val && val > 0.0 {
            max_mono_val = val;
            mono_nucleosomal_peak = Some(x);
        }
    }

    // Find di-nucleosomal peak (300-450 bp range)
    let mut di_nucleosomal_peak = None;
    let mut max_di_val = -1.0;
    for x in 300..=450 {
        let val = get_smoothed_count(x);
        if val > max_di_val && val > 0.0 {
            max_di_val = val;
            di_nucleosomal_peak = Some(x);
        }
    }

    let cfdna_ratio = if mono_cfdna_count > 0 {
        Some(short_cfdna_count as f64 / mono_cfdna_count as f64)
    } else {
        None
    };

    let mono_di_ratio = if di_nucleosomal_count > 0 {
        Some(mono_nucleosomal_count as f64 / di_nucleosomal_count as f64)
    } else {
        None
    };

    let short_fraction = sub_nucleosomal_count as f64 / total_inserts as f64;
    let mono_fraction = mono_nucleosomal_count as f64 / total_inserts as f64;
    let di_fraction = di_nucleosomal_count as f64 / total_inserts as f64;

    Some(InsertSizeMetrics {
        total_inserts,
        mono_nucleosomal_peak,
        di_nucleosomal_peak,
        mono_nucleosomal_count,
        di_nucleosomal_count,
        sub_nucleosomal_count,
        short_cfdna_count,
        mono_cfdna_count,
        cfdna_ratio,
        mono_di_ratio,
        short_fraction,
        mono_fraction,
        di_fraction,
        large_fragment_count,
    })
}

fn write_insert_size_metrics(config: &AlignmentQcConfig, metrics: &InsertSizeMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.align.insert_size_metrics.tsv", config.output_prefix))?;
    writeln!(out, "metric\tvalue")?;
    writeln!(out, "total_inserts\t{}", metrics.total_inserts)?;
    writeln!(
        out,
        "mono_nucleosomal_peak\t{}",
        metrics.mono_nucleosomal_peak.map(|x| x.to_string()).unwrap_or_else(|| "NA".to_string())
    )?;
    writeln!(
        out,
        "di_nucleosomal_peak\t{}",
        metrics.di_nucleosomal_peak.map(|x| x.to_string()).unwrap_or_else(|| "NA".to_string())
    )?;
    writeln!(out, "mono_nucleosomal_count\t{}", metrics.mono_nucleosomal_count)?;
    writeln!(out, "di_nucleosomal_count\t{}", metrics.di_nucleosomal_count)?;
    writeln!(out, "sub_nucleosomal_count\t{}", metrics.sub_nucleosomal_count)?;
    writeln!(out, "short_cfdna_count\t{}", metrics.short_cfdna_count)?;
    writeln!(out, "mono_cfdna_count\t{}", metrics.mono_cfdna_count)?;
    writeln!(
        out,
        "cfdna_ratio\t{}",
        metrics.cfdna_ratio.map(|v| format!("{v:.6}")).unwrap_or_else(|| "NA".to_string())
    )?;
    writeln!(
        out,
        "mono_di_ratio\t{}",
        metrics.mono_di_ratio.map(|v| format!("{v:.6}")).unwrap_or_else(|| "NA".to_string())
    )?;
    writeln!(out, "short_fraction\t{:.6}", metrics.short_fraction)?;
    writeln!(out, "mono_fraction\t{:.6}", metrics.mono_fraction)?;
    writeln!(out, "di_fraction\t{:.6}", metrics.di_fraction)?;
    writeln!(out, "large_fragment_count\t{}", metrics.large_fragment_count)?;
    Ok(())
}

pub fn sample_name_from_alignment_path(path: &str) -> String {
    Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(path)
        .trim_end_matches(".bam")
        .trim_end_matches(".cram")
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accuracy_quantiles_use_one_minus_de() {
        let mut acc = AccuracyMetrics::default();
        acc.observe(0.02);
        acc.observe(0.01);
        acc.observe(0.05);

        assert_eq!(acc.records_with_de, 3);
        assert!((acc.mean().unwrap() - 0.973333).abs() < 0.00001);
        assert_eq!(acc.median().unwrap(), 0.98);
    }

    #[test]
    fn accuracy_mode_uses_highest_histogram_bin_count() {
        let mut acc = AccuracyMetrics::default();
        acc.observe(0.01);
        acc.observe(0.01);
        acc.observe(0.03);

        assert_eq!(acc.mode().unwrap(), 0.99);
    }

    #[test]
    fn test_calculate_insert_size_metrics() {
        let mut hist = BTreeMap::new();
        // Mono-nucleosomal peak at 167 bp (natural bell-curve gradient)
        hist.insert(167, 100);
        hist.insert(166, 75);
        hist.insert(168, 75);
        hist.insert(165, 50);
        hist.insert(169, 50);
        hist.insert(164, 20);
        hist.insert(170, 20);

        // Di-nucleosomal peak at 330 bp (natural bell-curve gradient)
        hist.insert(330, 50);
        hist.insert(329, 30);
        hist.insert(331, 30);
        hist.insert(328, 15);
        hist.insert(332, 15);
        hist.insert(327, 5);
        hist.insert(333, 5);

        // Sub-nucleosomal / short cfDNA
        hist.insert(120, 20); // sub-nucleosome (30..143) and short cfDNA (100..150)
        hist.insert(80, 10);  // sub-nucleosome (30..143)
        hist.insert(1100, 5); // large (> 1000)

        let res_opt = calculate_insert_size_metrics(&hist);
        assert!(res_opt.is_some());
        let res = res_opt.unwrap();

        // Check peaks (smoothed)
        assert_eq!(res.mono_nucleosomal_peak, Some(167));
        assert_eq!(res.di_nucleosomal_peak, Some(330));

        // Check counts
        // mono-nucleosomal (143..=220): 100 + 75*2 + 50*2 + 20*2 = 390
        assert_eq!(res.mono_nucleosomal_count, 390);
        // di-nucleosomal (320..=480): 50 + 30*2 + 15*2 + 5*2 = 150
        assert_eq!(res.di_nucleosomal_count, 150);
        // sub-nucleosomal (30..143): 20 (at 120) + 10 (at 80) = 30
        assert_eq!(res.sub_nucleosomal_count, 30);
        // large (> 1000): 5
        assert_eq!(res.large_fragment_count, 5);

        // short cfDNA (100..150): 20
        assert_eq!(res.short_cfdna_count, 20);
        // mono cfDNA (151..220): 390
        assert_eq!(res.mono_cfdna_count, 390);

        // ratios
        assert!((res.cfdna_ratio.unwrap() - 20.0 / 390.0).abs() < 1e-6); // 20 / 390
        assert!((res.mono_di_ratio.unwrap() - 390.0 / 150.0).abs() < 1e-6);
    }
}
