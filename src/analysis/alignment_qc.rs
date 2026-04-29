use crate::analysis::bam_scan::{scan_bam_stream, BamScanConfig};
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
    fn observe_record(&mut self, header: &sam::Header, record: &bam::Record, mapq_threshold: u8) {
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
                let label = cigar_label(op.kind());
                *self.cigar_op_bases.entry(label.to_string()).or_insert(0) += len;
                match op.kind() {
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

        let contig = reference_name(header, record).unwrap_or_else(|| {
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

    scan_bam_stream(
        path,
        &BamScanConfig {
            threads: config.threads,
            show_progress: config.show_progress,
        },
        metrics,
        |state, header, record| {
            state.observe_record(header, record, config.mapq_threshold);
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

fn reference_name(header: &sam::Header, record: &bam::Record) -> Option<String> {
    let id = record.reference_sequence_id()?.ok()?;
    header
        .reference_sequences()
        .get_index(usize::from(id))
        .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).to_string())
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
}
