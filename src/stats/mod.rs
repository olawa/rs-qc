pub mod plotting;
use crate::analysis::index::AnnotationIndex;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::Ordering;

pub fn fraction(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

#[derive(Debug)]
pub struct AggregatedStats {
    pub dist_3p_means: Vec<f64>,
    pub dist_3p_sums_raw_all: Vec<f64>,
    pub dist_3p_sums_raw_analyzed: Vec<f64>,
    pub dist_3p_support: Vec<usize>,
    pub percentile_means: Vec<f64>,
    pub percentile_normalized: Vec<f64>,
    pub percentile_support: Vec<usize>,
    pub stratified_percentile_normalized: HashMap<String, Vec<f64>>,
    pub bin_size: usize,
    pub active_genes: usize,
    pub active_3p_genes: usize,
    pub total_reads: u64,
    pub gene_qc: Vec<ThreePrimeGeneQc>,
}

#[derive(Debug, Clone)]
pub struct ThreePrimeGeneQc {
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub transcript_id: Option<String>,
    pub chrom: String,
    pub strand: char,
    pub gene_len: u32,
    pub total_pct_count: u32,
    pub anchor_sum: u64,
    pub anchor_mean: f64,
    pub anchor_nonzero_bins: usize,
    pub max_ratio: f64,
    pub used: bool,
    pub skip_reason: Option<String>,
    pub counts_3p: Vec<u32>,
    pub counts_percentile: Vec<u32>,
}

pub struct ThreePrimeParams {
    pub normalization_bp: usize,
    pub max_3p_dist: usize,
    pub bin_size: usize,
    pub min_anchor_count: u64,
    pub min_anchor_mean: f64,
    pub min_anchor_nonzero_bins: usize,
    pub max_ratio: f64,
}

pub fn aggregate_genes(
    index: &AnnotationIndex,
    min_support: usize,
    three_prime_params: &ThreePrimeParams,
    total_reads: u64,
    _is_ends_mode: bool,
) -> AggregatedStats {
    let n_3p_bins = (three_prime_params.max_3p_dist + three_prime_params.bin_size - 1)
        / three_prime_params.bin_size;
    let norm_bins = (three_prime_params.normalization_bp + three_prime_params.bin_size - 1)
        / three_prime_params.bin_size;

    let (
        dist_3p_sums,
        dist_3p_support,
        dist_3p_sums_raw,
        percentile_sums,
        percentile_support,
        active_genes,
        active_3p_genes,
        dist_3p_raw_analyzed,
        gene_qc,
    ) = index
        .genes
        .par_iter()
        .fold(
            || {
                (
                    vec![0.0f64; n_3p_bins],
                    vec![0usize; n_3p_bins],
                    vec![0.0f64; n_3p_bins],
                    vec![0.0f64; 100],
                    vec![0usize; 100],
                    0usize,
                    0usize,
                    vec![0.0f64; n_3p_bins],
                    Vec::new(),
                )
            },
            |(
                mut d3p_s,
                mut d3p_sup,
                mut d3p_raw,
                mut p_s,
                mut p_sup,
                mut active,
                mut active_3p,
                mut d3p_raw_a,
                mut qc_list,
            ),
             gene| {
                let total_pct_count: u32 = gene
                    .counts_percentile
                    .iter()
                    .map(|a| a.load(Ordering::Relaxed))
                    .sum();

                if total_pct_count >= min_support as u32 {
                    active += 1;
                    let mean_pct = (total_pct_count as f64 / 100.0).max(0.001);
                    for (i, count_atomic) in gene.counts_percentile.iter().enumerate() {
                        let count = count_atomic.load(Ordering::Relaxed);
                        p_s[i] += count as f64 / mean_pct;
                        p_sup[i] += 1;
                    }
                }

                let gene_len = gene.total_len;
                let valid_3p_bins = gene.counts_3p.len().min(n_3p_bins);
                let anchor_bins = norm_bins.min(valid_3p_bins).max(1);

                let mut anchor_sum = 0u64;
                let mut anchor_nonzero = 0usize;
                for i in 0..anchor_bins {
                    let c = gene.counts_3p[i].load(Ordering::Relaxed) as u64;
                    anchor_sum += c;
                    if c > 0 {
                        anchor_nonzero += 1;
                    }
                }

                let anchor_mean = anchor_sum as f64 / anchor_bins as f64;

                let mut used = true;
                let mut skip_reason = None;

                if total_pct_count < min_support as u32 {
                    used = false;
                    skip_reason = Some(format!("low total support ({})", total_pct_count));
                } else if gene_len < three_prime_params.normalization_bp as u32 {
                    used = false;
                    skip_reason = Some(format!("too short ({})", gene_len));
                } else if anchor_sum < three_prime_params.min_anchor_count {
                    used = false;
                    skip_reason = Some(format!("low anchor sum ({})", anchor_sum));
                } else if anchor_mean < three_prime_params.min_anchor_mean {
                    used = false;
                    skip_reason = Some(format!("low anchor mean ({:.2})", anchor_mean));
                } else if anchor_nonzero < three_prime_params.min_anchor_nonzero_bins {
                    used = false;
                    skip_reason = Some(format!("low anchor nonzero bins ({})", anchor_nonzero));
                }

                let mut max_ratio_found: f64 = 0.0;
                if used {
                    active_3p += 1;
                    for i in 0..valid_3p_bins {
                        if (i * three_prime_params.bin_size) as u32 >= gene_len {
                            break;
                        }
                        let count = gene.counts_3p[i].load(Ordering::Relaxed);
                        let ratio = count as f64 / anchor_mean;
                        max_ratio_found = max_ratio_found.max(ratio);

                        let clipped_ratio = ratio.min(three_prime_params.max_ratio);
                        d3p_s[i] += clipped_ratio;
                        d3p_sup[i] += 1;
                        d3p_raw_a[i] += count as f64;
                    }
                }

                // Add to raw all plot
                for i in 0..valid_3p_bins {
                    if (i * three_prime_params.bin_size) as u32 >= gene_len {
                        break;
                    }
                    d3p_raw[i] += gene.counts_3p[i].load(Ordering::Relaxed) as f64;
                }

                qc_list.push(ThreePrimeGeneQc {
                    gene_id: gene.id.clone(),
                    gene_name: gene.name.clone(),
                    transcript_id: Some(gene.representative.id.clone()),
                    chrom: gene.chrom.clone(),
                    strand: gene.representative.strand,
                    gene_len,
                    total_pct_count,
                    anchor_sum,
                    anchor_mean,
                    anchor_nonzero_bins: anchor_nonzero,
                    max_ratio: max_ratio_found,
                    used,
                    skip_reason,
                    counts_3p: gene
                        .counts_3p
                        .iter()
                        .map(|a| a.load(Ordering::Relaxed))
                        .collect(),
                    counts_percentile: gene
                        .counts_percentile
                        .iter()
                        .map(|a| a.load(Ordering::Relaxed))
                        .collect(),
                });

                (
                    d3p_s, d3p_sup, d3p_raw, p_s, p_sup, active, active_3p, d3p_raw_a, qc_list,
                )
            },
        )
        .reduce(
            || {
                (
                    vec![0.0f64; n_3p_bins],
                    vec![0usize; n_3p_bins],
                    vec![0.0f64; n_3p_bins],
                    vec![0.0f64; 100],
                    vec![0usize; 100],
                    0usize,
                    0usize,
                    vec![0.0f64; n_3p_bins],
                    Vec::new(),
                )
            },
            |mut a, mut b| {
                for i in 0..n_3p_bins {
                    a.0[i] += b.0[i];
                    a.1[i] += b.1[i];
                    a.2[i] += b.2[i];
                    a.7[i] += b.7[i];
                }
                for i in 0..100 {
                    a.3[i] += b.3[i];
                    a.4[i] += b.4[i];
                }
                a.5 += b.5;
                a.6 += b.6;
                a.8.append(&mut b.8);
                a
            },
        );

    let dist_3p_means: Vec<f64> = dist_3p_sums
        .into_iter()
        .enumerate()
        .map(|(i, sum)| {
            if dist_3p_support[i] > 0 {
                sum / dist_3p_support[i] as f64
            } else {
                0.0
            }
        })
        .collect();

    let percentile_means: Vec<f64> = percentile_sums
        .into_iter()
        .enumerate()
        .map(|(i, s)| s / percentile_support[i].max(1) as f64)
        .collect();
    let max_percentile = percentile_means.iter().cloned().fold(0.0_f64, f64::max);
    let percentile_normalized = if max_percentile > 0.0 {
        percentile_means
            .iter()
            .map(|&v| (v / max_percentile) * 100.0)
            .collect()
    } else {
        vec![0.0; 100]
    };

    AggregatedStats {
        dist_3p_means,
        dist_3p_sums_raw_all: dist_3p_sums_raw,
        dist_3p_sums_raw_analyzed: dist_3p_raw_analyzed,
        dist_3p_support,
        percentile_means,
        percentile_normalized,
        percentile_support,
        stratified_percentile_normalized: HashMap::new(),
        bin_size: three_prime_params.bin_size,
        active_genes,
        active_3p_genes,
        total_reads,
        gene_qc,
    }
}

/// Computes gene body coverage using the RSeQC-compatible approach:
///
/// **Default (RSeQC-compatible, `coverage_weighted = false`):**
/// 1. Filter genes with `total_bin_count < min_support`.
/// 2. Per-gene normalize: `profile[i] = count[i] / mean_depth` (mean bin depth → 1.0).
///    This ensures a high-coverage short transcript contributes its *shape*, not its
///    *amplitude*, so it never dominates the average over lowly-expressed long genes.
/// 3. Simple equal-weight mean across all qualifying genes.
/// 4. Scale to 0–100 % relative to peak.
///
/// **Optional weighted mode (`coverage_weighted = true`):**
/// Uses `w = min(total_count, 20 × min_support)` as the per-gene weight, so
/// well-covered genes have more influence than barely-expressed ones.
/// Useful for exploratory analysis; not recommended as a primary QC metric.
///
/// Returns `(raw_means, normalized_0_to_100)`.
pub fn aggregate_rseqc_classic(
    index: &AnnotationIndex,
    _is_ends_mode: bool,
    min_support: u32,
    coverage_weighted: bool,
) -> (Vec<f64>, Vec<f64>) {
    let weight_cap = (min_support as f64 * 20.0).max(2000.0);

    let (profile_sums, total_weight, total_genes) = index
        .genes
        .par_iter()
        .fold(
            || (vec![0.0f64; 100], 0.0f64, 0usize),
            |(mut sums, mut tw, mut n), gene| {
                if gene.total_len < 100 {
                    return (sums, tw, n);
                }
                let total_count: u32 = gene
                    .counts_percentile
                    .iter()
                    .map(|a| a.load(Ordering::Relaxed))
                    .sum();
                if total_count < min_support as u32 {
                    return (sums, tw, n);
                }
                let mean_depth = (total_count as f64) / 100.0;
                let weight = if coverage_weighted {
                    (total_count as f64).min(weight_cap)
                } else {
                    1.0
                };
                for (i, count_atomic) in gene.counts_percentile.iter().enumerate() {
                    let count = count_atomic.load(Ordering::Relaxed);
                    sums[i] += weight * (count as f64 / mean_depth);
                }
                tw += weight;
                n += 1;
                (sums, tw, n)
            },
        )
        .reduce(
            || (vec![0.0f64; 100], 0.0f64, 0usize),
            |mut a, b| {
                for i in 0..100 {
                    a.0[i] += b.0[i];
                }
                a.1 += b.1;
                a.2 += b.2;
                a
            },
        );

    if total_genes == 0 || total_weight == 0.0 {
        return (vec![0.0; 100], vec![0.0; 100]);
    }

    println!(
        "  - Gene body coverage: {} genes passed min_support={}{}",
        total_genes,
        min_support,
        if coverage_weighted {
            " (coverage-weighted)"
        } else {
            ""
        }
    );

    let raw_means: Vec<f64> = profile_sums.iter().map(|&s| s / total_weight).collect();

    // Scale to 0–100 % by peak value.
    let mut normalized = raw_means.clone();
    let max_val = normalized.iter().cloned().fold(0.0_f64, f64::max);
    if max_val > 0.0 {
        for val in normalized.iter_mut() {
            *val = (*val / max_val) * 100.0;
        }
    }

    (raw_means, normalized)
}

/// Length-stratified gene body coverage.
///
/// Applies the same per-gene normalization as [`aggregate_rseqc_classic`] independently
/// within each transcript-length class. Returns a `HashMap` whose keys are the class
/// labels and whose values are 100-element vectors normalized to 0–100 % by peak.
///
/// The returned map feeds directly into [`plotting::generate_gene_body_plot`], producing
/// a single multi-series SVG that makes it easy to spot whether 5' bias is concentrated
/// in a particular length range (a hallmark of long-3'-UTR artefacts or degradation).
///
/// Length classes (spliced transcript length, in bp):
///
/// | Class          | Range        | Notes                                           |
/// |----------------|--------------|-------------------------------------------------|
/// | `short`        | < 1,500      | Compact housekeeping genes; minimal UTR issues  |
/// | `medium`       | 1,500–5,000  | Typical protein-coding transcripts              |
/// | `long`         | > 5,000      | Long UTRs; most susceptible to 5' artefacts     |
pub fn aggregate_rseqc_stratified(
    index: &AnnotationIndex,
    min_support: u32,
    coverage_weighted: bool,
) -> HashMap<String, Vec<f64>> {
    // (label, lower_inclusive_bp, upper_exclusive_bp)
    const BINS: &[(&str, u32, u32)] = &[
        ("short (<1.5kb)", 0, 1_500),
        ("medium (1.5-5kb)", 1_500, 5_001),
        ("long (5-10kb)", 5_001, 10_001),
        ("very long (>10kb)", 10_001, u32::MAX),
    ];

    let weight_cap = (min_support as f64 * 20.0).max(2000.0);

    // One accumulator per bin: (profile_sums[100], total_weight)
    let accum: Vec<(Vec<f64>, f64, usize)> = index
        .genes
        .par_iter()
        .fold(
            || vec![(vec![0.0f64; 100], 0.0f64, 0usize); BINS.len()],
            |mut local_accum, gene| {
                let tx_len = gene.total_len;
                if let Some(bin_idx) = BINS
                    .iter()
                    .position(|&(_, lo, hi)| tx_len >= lo && tx_len < hi)
                {
                    let total_count: u32 = gene
                        .counts_percentile
                        .iter()
                        .map(|a| a.load(Ordering::Relaxed))
                        .sum();
                    if total_count >= min_support as u32 {
                        let mean_depth = (total_count as f64) / 100.0;
                        let weight = if coverage_weighted {
                            (total_count as f64).min(weight_cap)
                        } else {
                            1.0
                        };
                        let (sums, tw, n) = &mut local_accum[bin_idx];
                        for (i, count_atomic) in gene.counts_percentile.iter().enumerate() {
                            let count = count_atomic.load(Ordering::Relaxed);
                            sums[i] += weight * (count as f64 / mean_depth);
                        }
                        *tw += weight;
                        *n += 1;
                    }
                }
                local_accum
            },
        )
        .reduce(
            || vec![(vec![0.0f64; 100], 0.0f64, 0usize); BINS.len()],
            |mut a, b| {
                for i in 0..a.len() {
                    for j in 0..100 {
                        a[i].0[j] += b[i].0[j];
                    }
                    a[i].1 += b[i].1;
                    a[i].2 += b[i].2;
                }
                a
            },
        );

    let mut result = HashMap::new();
    for (i, &(label, _, _)) in BINS.iter().enumerate() {
        let (sums, total_weight, n_genes) = &accum[i];
        if *total_weight == 0.0 {
            continue;
        }
        println!(
            "  - {:20} {:>5} genes (min_support={})",
            label, n_genes, min_support
        );
        let raw: Vec<f64> = sums.iter().map(|&s| s / total_weight).collect();
        let max_val = raw.iter().cloned().fold(0.0_f64, f64::max);
        let normalized: Vec<f64> = if max_val > 0.0 {
            raw.iter().map(|&v| (v / max_val) * 100.0).collect()
        } else {
            vec![0.0; 100]
        };
        result.insert(label.to_string(), normalized);
    }
    result
}

pub fn write_tsv(path: &str, header: &str, data: &[f64], support: &[usize]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "{}", header)?;
    for i in 0..data.len() {
        writeln!(file, "{}\t{:.4}\t{}", i, data[i], support[i])?;
    }
    Ok(())
}

pub fn write_3p_wide_tsv<T: std::fmt::Display>(
    path: &str,
    results: &HashMap<String, Vec<T>>,
    bin_size: usize,
    max_3p_dist: usize,
) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    let mut sample_ids: Vec<_> = results.keys().collect();
    sample_ids.sort();

    write!(file, "distance_start_bp\tdistance_end_bp")?;
    for sample_id in &sample_ids {
        write!(file, "\t{}", sample_id)?;
    }
    writeln!(file)?;

    if sample_ids.is_empty() {
        return Ok(());
    }

    let n_bins = results
        .get(sample_ids[0])
        .map(|values| values.len())
        .unwrap_or(0);

    for i in 0..n_bins {
        let start = i * bin_size;
        let end = ((i + 1) * bin_size).min(max_3p_dist);
        write!(file, "{}\t{}", start, end)?;
        for sample_id in &sample_ids {
            let values = results
                .get(*sample_id)
                .expect("sample missing from results");
            write!(file, "\t{}", values[i])?;
        }
        writeln!(file)?;
    }

    Ok(())
}

pub fn write_classic_wide_format(
    path: &str,
    results: &HashMap<String, Vec<f64>>,
) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    // MultiQC header: Percentile	1	2	3	...	100
    write!(file, "Percentile")?;
    for i in 1..=100 {
        write!(file, "\t{}", i)?;
    }
    writeln!(file)?;

    for (sample_id, data) in results {
        write!(file, "{}", sample_id)?;
        for val in data {
            write!(file, "\t{:.4}", val)?;
        }
        writeln!(file)?;
    }
    Ok(())
}

pub fn write_gene_profiles_tsv(path: &str, qc_data: &[ThreePrimeGeneQc]) -> anyhow::Result<()> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "gene_id\tgene_name\ttranscript_id\tchrom\tstrand\tlength\tgene_body_bin_count\tanchor_3p_count\tmean_body_coverage\tanchor_mean_coverage\tused_for_3p\tskip_reason\tcounts_3p_csv\tcounts_percentile_csv"
    )?;

    for q in qc_data {
        let counts_3p_csv = q
            .counts_3p
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let counts_pct_csv = q
            .counts_percentile
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(",");

        let mean_body = q.total_pct_count as f64 / 100.0;

        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
            q.gene_id,
            q.gene_name.as_deref().unwrap_or("-"),
            q.transcript_id.as_deref().unwrap_or("-"),
            q.chrom,
            q.strand,
            q.gene_len,
            q.total_pct_count,
            q.anchor_sum,
            mean_body,
            q.anchor_mean,
            q.used,
            q.skip_reason.as_deref().unwrap_or("-"),
            counts_3p_csv,
            counts_pct_csv
        )?;
    }

    Ok(())
}

pub fn write_feature_counts_tsv(path: &str, qc_data: &[ThreePrimeGeneQc]) -> anyhow::Result<()> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "gene_id\tgene_name\tchrom\tstrand\tlength\tgene_body_coverage_count\tanchor_3p_count"
    )?;

    for q in qc_data {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            q.gene_id,
            q.gene_name.as_deref().unwrap_or("-"),
            q.chrom,
            q.strand,
            q.gene_len,
            q.total_pct_count,
            q.anchor_sum
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_rseqc_sampling_indices() {
        // Verification against RSeQC formula: int(round(i * (L - 1.0) / 99.0))
        let lengths = vec![100.0, 1000.0, 500.0];
        for len in lengths {
            let mut indices = Vec::new();
            for i in 0..100 {
                let idx = (i as f64 * (len - 1.0) / 99.0).round() as usize;
                indices.push(idx);
            }
            assert_eq!(indices.len(), 100);
            assert_eq!(indices[0], 0);
            assert_eq!(indices[99], (len - 1.0) as usize);

            // Check specific midpoints
            if len == 100.0 {
                assert_eq!(indices[50], 50);
            }
            if len == 500.0 {
                // int(round(50 * 499 / 99)) = int(round(25202 / 99)) = int(round(252.04)) = 252
                assert_eq!(indices[50], 252);
            }
        }
    }

    #[test]
    fn test_3p_robustness_and_normalization() {
        use super::*;
        use crate::analysis::index::AnnotationIndex;
        use crate::models::transcript::{Exon, Transcript};
        use crate::models::Gene;
        use std::sync::atomic::Ordering;

        let params = ThreePrimeParams {
            normalization_bp: 200,
            max_3p_dist: 1000,
            bin_size: 100,
            min_anchor_count: 50,
            min_anchor_mean: 5.0,
            min_anchor_nonzero_bins: 2,
            max_ratio: 3.0,
        };

        // 1. Accepted Gene: anchor bins 0, 1 (200bp total).
        // Counts: bin0=100, bin1=100. Anchor mean = 100.
        let tx1 = Transcript::new(
            "tx1".into(),
            "chr1".into(),
            '+',
            None,
            vec![Exon {
                start: 1,
                end: 1001,
            }],
            None,
            None,
        );
        let gene1 = Gene::new("g1".into(), None, tx1, 1000, 100);
        gene1.counts_3p[0].store(100, Ordering::Relaxed);
        gene1.counts_3p[1].store(100, Ordering::Relaxed);
        gene1.counts_3p[2].store(50, Ordering::Relaxed);
        for i in 0..100 {
            gene1.counts_percentile[i].store(100, Ordering::Relaxed);
        }

        // 2. Low Anchor Gene: bin0=1, bin1=1 (anchor sum=2 < 50). Should be skipped.
        let tx2 = Transcript::new(
            "tx2".into(),
            "chr1".into(),
            '+',
            None,
            vec![Exon {
                start: 1,
                end: 1001,
            }],
            None,
            None,
        );
        let gene2 = Gene::new("g2".into(), None, tx2, 1000, 100);
        gene2.counts_3p[0].store(1, Ordering::Relaxed);
        gene2.counts_3p[1].store(1, Ordering::Relaxed);
        gene2.counts_3p[2].store(1000, Ordering::Relaxed); // Huge spike downstream
        for i in 0..100 {
            gene2.counts_percentile[i].store(100, Ordering::Relaxed);
        }

        // 3. High Ratio Gene: anchor sum=200, mean=100. bin2=500 -> ratio 5.0. Should be clipped to 3.0.
        let tx3 = Transcript::new(
            "tx3".into(),
            "chr1".into(),
            '+',
            None,
            vec![Exon {
                start: 1,
                end: 1001,
            }],
            None,
            None,
        );
        let gene3 = Gene::new("g3".into(), None, tx3, 1000, 100);
        gene3.counts_3p[0].store(100, Ordering::Relaxed);
        gene3.counts_3p[1].store(100, Ordering::Relaxed);
        gene3.counts_3p[2].store(500, Ordering::Relaxed); // Ratio 5.0
        for i in 0..100 {
            gene3.counts_percentile[i].store(100, Ordering::Relaxed);
        }

        let index = AnnotationIndex {
            version: 1,
            genes: std::sync::Arc::new(vec![gene1, gene2, gene3]),
            feature_sizes: Default::default(),
            chrom_spans: Default::default(),
            genes_by_chrom: Default::default(),
        };

        let stats = aggregate_genes(&index, 10, &params, 10000, false);

        // Gene 1 and 3 are kept. Gene 2 is skipped.
        assert_eq!(stats.active_3p_genes, 2);

        // Gene 1: bin0=1.0, bin1=1.0, bin2=0.5
        // Gene 3: bin0=1.0, bin1=1.0, bin2=3.0 (clipped from 5.0)
        // Average: bin0=1.0, bin1=1.0, bin2=(0.5 + 3.0)/2 = 1.75

        assert!((stats.dist_3p_means[0] - 1.0).abs() < 0.001);
        assert!((stats.dist_3p_means[1] - 1.0).abs() < 0.001);
        assert!((stats.dist_3p_means[2] - 1.75).abs() < 0.001);

        // Check QC data
        let g2_qc = stats.gene_qc.iter().find(|q| q.gene_id == "g2").unwrap();
        assert!(!g2_qc.used);
        assert!(g2_qc
            .skip_reason
            .as_ref()
            .unwrap()
            .contains("low anchor sum"));
    }
}
