use crate::io::bam::{reference_span, ReferenceNames};
use crate::io::text::open_maybe_gz;
use crate::models::normalize_chrom;
use crate::stats::fraction;
use anyhow::{Context, Result};
use noodles::bam;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ContaminationQcConfig {
    pub inputs: Vec<String>,
    pub output_prefix: String,
    pub rdna_bed: Option<String>,
    pub rdna_contigs: Vec<String>,
    pub rrna_fasta: Option<String>,
    pub kmer_size: usize,
    pub min_kmer_hits: usize,
    pub sample_size: usize,
    pub mapq_threshold: u8,
    pub scan_all_reads: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct ContaminationQcSummary {
    pub sample: String,
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub mtdna_reads: u64,
    pub rdna_interval_reads: u64,
    pub rdna_contig_reads: u64,
    pub rrna_kmer_reads: u64,
    pub rrna_kmer_reads_screened: u64,
    pub mtdna_fraction: f64,
    pub rdna_aligned_fraction: f64,
    pub rrna_kmer_fraction: f64,
}

#[derive(Debug, Default)]
struct ContaminationCounts {
    total_reads: u64,
    mapped_reads: u64,
    mtdna_reads: u64,
    rdna_interval_reads: u64,
    rdna_contig_reads: u64,
    rdna_aligned_reads: u64,
    rrna_kmer_reads: u64,
    rrna_kmer_reads_screened: u64,
}

pub fn run_contamination_qc(config: &ContaminationQcConfig) -> Result<Vec<ContaminationQcSummary>> {
    let rdna_intervals = if let Some(path) = &config.rdna_bed {
        Some(ContaminantIndex::from_bed(path)?)
    } else {
        None
    };
    let rdna_contigs = normalized_name_set(&config.rdna_contigs);
    let rrna_kmers = if let Some(path) = &config.rrna_fasta {
        Some(build_kmer_index(path, config.kmer_size)?)
    } else {
        None
    };

    let mut summaries = Vec::new();
    for input in &config.inputs {
        let summary = scan_bam_contamination(
            input,
            rdna_intervals.as_ref(),
            &rdna_contigs,
            rrna_kmers.as_ref(),
            config,
        )?;
        write_contamination_outputs(config, &summary)?;
        summaries.push(summary);
    }
    Ok(summaries)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContaminantInterval {
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromContaminantIndex {
    pub intervals: Arc<[ContaminantInterval]>,
}

#[derive(Debug, Default, Clone)]
pub struct ContaminantCursor {
    idx: usize,
}

#[derive(Debug, Default, Clone)]
pub struct ContaminantIndex {
    pub chroms: HashMap<String, Arc<ChromContaminantIndex>>,
}

impl ContaminantIndex {
    pub fn from_bed(path: &str) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("Failed to open rDNA BED file: {}", path))?;
        let reader = BufReader::new(file);
        let mut by_chrom: HashMap<String, Vec<ContaminantInterval>> = HashMap::new();

        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                anyhow::bail!("Invalid BED line {} in {}", line_no + 1, path);
            }

            let chrom = normalize_chrom(fields[0]);
            let start: u64 = fields[1]
                .parse()
                .with_context(|| format!("Invalid BED start at line {}", line_no + 1))?;
            let end: u64 = fields[2]
                .parse()
                .with_context(|| format!("Invalid BED end at line {}", line_no + 1))?;

            if start < end {
                by_chrom
                    .entry(chrom.into_owned())
                    .or_default()
                    .push(ContaminantInterval { start, end });
            }
        }

        let chroms = by_chrom
            .into_iter()
            .filter_map(|(chrom, intervals)| {
                let merged = merge_intervals(intervals);
                if merged.is_empty() {
                    None
                } else {
                    Some((
                        chrom,
                        Arc::new(ChromContaminantIndex {
                            intervals: merged.into(),
                        }),
                    ))
                }
            })
            .collect();

        Ok(Self { chroms })
    }
}

impl ChromContaminantIndex {
    pub fn cursor_at(&self, pos: u64) -> ContaminantCursor {
        ContaminantCursor {
            idx: self.intervals.partition_point(|iv| iv.end <= pos),
        }
    }

    pub fn contains(&self, pos: u64, cursor: &mut ContaminantCursor) -> bool {
        while cursor.idx < self.intervals.len() && self.intervals[cursor.idx].end <= pos {
            cursor.idx += 1;
        }

        self.intervals
            .get(cursor.idx)
            .is_some_and(|iv| iv.start <= pos && pos < iv.end)
    }
}

fn scan_bam_contamination(
    path: &str,
    rdna_intervals: Option<&ContaminantIndex>,
    rdna_contigs: &HashSet<String>,
    rrna_kmers: Option<&HashSet<Vec<u8>>>,
    config: &ContaminationQcConfig,
) -> Result<ContaminationQcSummary> {
    let file = File::open(path).with_context(|| format!("failed to open BAM {path}"))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header()?;
    let ref_names = ReferenceNames::new(&header);
    let mut counts = ContaminationCounts::default();

    for result in reader.records() {
        let record = result?;
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        counts.total_reads += 1;
        if !flags.is_unmapped() {
            counts.mapped_reads += 1;
        }

        if let Some(norm) = ref_names.normalized(&record) {
            if is_mtdna_chrom(norm) {
                counts.mtdna_reads += 1;
            }
            let mut is_rdna_aligned = false;
            if rdna_contigs.contains(norm) {
                counts.rdna_contig_reads += 1;
                is_rdna_aligned = true;
            }
            if let Some(index) = rdna_intervals {
                if read_overlaps_contaminant_interval(index, norm, &record) {
                    counts.rdna_interval_reads += 1;
                    is_rdna_aligned = true;
                }
            }
            if is_rdna_aligned {
                counts.rdna_aligned_reads += 1;
            }
        }

        if let Some(kmers) = rrna_kmers {
            if should_kmer_screen(&record, config)
                && (config.sample_size == 0
                    || counts.rrna_kmer_reads_screened < config.sample_size as u64)
            {
                counts.rrna_kmer_reads_screened += 1;
                if sequence_has_kmer_hits_bam(
                    &record.sequence(),
                    kmers,
                    config.kmer_size,
                    config.min_kmer_hits,
                ) {
                    counts.rrna_kmer_reads += 1;
                }
            }
        }
    }

    let sample = std::path::Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(path)
        .to_string();
    Ok(summarize_contamination_counts(sample, counts))
}

fn summarize_contamination_counts(
    sample: String,
    counts: ContaminationCounts,
) -> ContaminationQcSummary {
    ContaminationQcSummary {
        sample,
        total_reads: counts.total_reads,
        mapped_reads: counts.mapped_reads,
        mtdna_reads: counts.mtdna_reads,
        rdna_interval_reads: counts.rdna_interval_reads,
        rdna_contig_reads: counts.rdna_contig_reads,
        rrna_kmer_reads: counts.rrna_kmer_reads,
        rrna_kmer_reads_screened: counts.rrna_kmer_reads_screened,
        mtdna_fraction: fraction(counts.mtdna_reads, counts.total_reads),
        rdna_aligned_fraction: fraction(counts.rdna_aligned_reads, counts.total_reads),
        rrna_kmer_fraction: fraction(counts.rrna_kmer_reads, counts.rrna_kmer_reads_screened),
    }
}

fn read_overlaps_contaminant_interval(
    index: &ContaminantIndex,
    chrom: &str,
    record: &bam::Record,
) -> bool {
    let Some(chrom_index) = index.chroms.get(chrom) else {
        return false;
    };
    let Some((start, end)) = reference_span(record) else {
        return false;
    };
    let idx = chrom_index.intervals.partition_point(|iv| iv.end <= start);
    chrom_index
        .intervals
        .get(idx)
        .is_some_and(|iv| iv.start < end && start < iv.end)
}

fn should_kmer_screen(record: &bam::Record, config: &ContaminationQcConfig) -> bool {
    if config.scan_all_reads {
        return true;
    }
    if record.flags().is_unmapped() {
        return true;
    }
    record
        .mapping_quality()
        .map(|mapq| mapq.get() < config.mapq_threshold)
        .unwrap_or(true)
}

fn normalized_name_set(values: &[String]) -> HashSet<String> {
    values
        .iter()
        .map(|v| normalize_chrom(v).into_owned())
        .collect()
}

fn is_mtdna_chrom(chrom: &str) -> bool {
    matches!(chrom, "m" | "mt" | "mitochondria")
}

fn write_contamination_outputs(
    config: &ContaminationQcConfig,
    summary: &ContaminationQcSummary,
) -> Result<()> {
    let path = format!(
        "{}.{}.contam.summary.tsv",
        config.output_prefix, summary.sample
    );
    let mut out = File::create(&path)?;
    writeln!(out, "metric\tvalue")?;
    writeln!(out, "total_reads\t{}", summary.total_reads)?;
    writeln!(out, "mapped_reads\t{}", summary.mapped_reads)?;
    writeln!(out, "mtdna_reads\t{}", summary.mtdna_reads)?;
    writeln!(out, "mtdna_fraction\t{:.6}", summary.mtdna_fraction)?;
    writeln!(out, "rdna_interval_reads\t{}", summary.rdna_interval_reads)?;
    writeln!(out, "rdna_contig_reads\t{}", summary.rdna_contig_reads)?;
    writeln!(
        out,
        "rdna_aligned_fraction\t{:.6}",
        summary.rdna_aligned_fraction
    )?;
    writeln!(out, "rrna_kmer_reads\t{}", summary.rrna_kmer_reads)?;
    writeln!(
        out,
        "rrna_kmer_reads_screened\t{}",
        summary.rrna_kmer_reads_screened
    )?;
    writeln!(out, "rrna_kmer_fraction\t{:.6}", summary.rrna_kmer_fraction)?;

    let json_path = format!(
        "{}.{}.contam.summary.json",
        config.output_prefix, summary.sample
    );
    let json = serde_json::to_string_pretty(summary)?;
    std::fs::write(json_path, json)?;
    Ok(())
}

fn build_kmer_index(path: &str, k: usize) -> Result<HashSet<Vec<u8>>> {
    anyhow::ensure!(k > 0, "--kmer-size must be greater than zero");
    let mut kmers = HashSet::new();
    let mut seq = Vec::new();
    for line in open_maybe_gz(path)?.lines() {
        let line = line?;
        if line.starts_with('>') {
            add_sequence_kmers(&seq, k, &mut kmers);
            seq.clear();
        } else {
            seq.extend(
                line.trim()
                    .as_bytes()
                    .iter()
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    add_sequence_kmers(&seq, k, &mut kmers);
    Ok(kmers)
}

fn add_sequence_kmers(seq: &[u8], k: usize, kmers: &mut HashSet<Vec<u8>>) {
    if seq.len() < k {
        return;
    }
    for window in seq.windows(k) {
        if window
            .iter()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        {
            kmers.insert(window.to_vec());
            kmers.insert(reverse_complement(window));
        }
    }
}

fn sequence_has_kmer_hits_bam(
    seq: &bam::record::Sequence,
    kmers: &HashSet<Vec<u8>>,
    k: usize,
    min_hits: usize,
) -> bool {
    if seq.len() < k || min_hits == 0 {
        return false;
    }
    let mut hits = 0_usize;
    let mut window = vec![0u8; k];
    let seq_bytes: Vec<u8> = seq.iter().map(u8::from).collect();
    for i in 0..=seq_bytes.len() - k {
        let mut valid = true;
        for j in 0..k {
            let b = seq_bytes[i + j].to_ascii_uppercase();
            if !matches!(b, b'A' | b'C' | b'G' | b'T') {
                valid = false;
                break;
            }
            window[j] = b;
        }
        if valid && kmers.contains(&window) {
            hits += 1;
            if hits >= min_hits {
                return true;
            }
        }
    }
    false
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            _ => b'N',
        })
        .collect()
}

fn merge_intervals(mut intervals: Vec<ContaminantInterval>) -> Vec<ContaminantInterval> {
    intervals.sort_by_key(|iv| (iv.start, iv.end));
    let mut merged: Vec<ContaminantInterval> = Vec::new();

    for interval in intervals {
        if let Some(last) = merged.last_mut() {
            if interval.start <= last.end {
                last.end = last.end.max(interval.end);
                continue;
            }
        }
        merged.push(interval);
    }

    merged
}

fn sequence_has_kmer_hits(
    seq: &[u8],
    kmers: &std::collections::HashSet<Vec<u8>>,
    k: usize,
    min_hits: usize,
) -> bool {
    if seq.len() < k || min_hits == 0 {
        return false;
    }
    let mut hits = 0;
    for window in seq.windows(k) {
        if kmers.contains(window) {
            hits += 1;
            if hits >= min_hits {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_intervals_collapses_overlaps() {
        let merged = merge_intervals(vec![
            ContaminantInterval { start: 20, end: 30 },
            ContaminantInterval { start: 10, end: 15 },
            ContaminantInterval { start: 14, end: 25 },
        ]);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 10);
        assert_eq!(merged[0].end, 30);
    }

    #[test]
    fn exact_kmer_screen_counts_reverse_complements() {
        let mut kmers = HashSet::new();
        add_sequence_kmers(b"ACGTAAAA", 4, &mut kmers);

        assert!(sequence_has_kmer_hits(b"TTTT", &kmers, 4, 1));
        assert!(sequence_has_kmer_hits(b"ACGT", &kmers, 4, 1));
        assert!(!sequence_has_kmer_hits(b"CCCC", &kmers, 4, 1));
    }

    #[test]
    fn rdna_aligned_fraction_uses_union_not_maximum_subset() {
        let counts = ContaminationCounts {
            total_reads: 10,
            rdna_interval_reads: 3,
            rdna_contig_reads: 4,
            rdna_aligned_reads: 7,
            ..ContaminationCounts::default()
        };

        let summary = summarize_contamination_counts("sample".to_string(), counts);

        assert_eq!(summary.rdna_interval_reads, 3);
        assert_eq!(summary.rdna_contig_reads, 4);
        assert!((summary.rdna_aligned_fraction - 0.7).abs() < 1e-9);
    }
}
