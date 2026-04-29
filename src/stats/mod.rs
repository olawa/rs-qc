pub mod plotting;
use crate::analysis::index::AnnotationIndex;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

pub struct AggregatedStats {
    pub dist_3p_means: Vec<f64>,
    pub dist_3p_support: Vec<usize>,
    pub percentile_means: Vec<f64>,
    pub percentile_support: Vec<usize>,
}

pub fn aggregate_genes(
    index: &AnnotationIndex,
    _normalization_bp: usize,
    min_support: usize,
    max_3p_dist: usize,
    is_ends_mode: bool,
) -> AggregatedStats {
    let mut dist_3p_sums = vec![0.0f64; max_3p_dist];
    let mut dist_3p_support = vec![0usize; max_3p_dist];

    let mut percentile_sums = vec![0.0f64; 100];
    let mut percentile_support = vec![0usize; 100];

    let mut active_genes = 0usize;

    // Explicitly iterate over the Arc<Vec<Gene>>
    for gene in index.genes.iter() {
        let tx_len = gene.total_len as u64;

        let mut coverage_3p = vec![0u32; max_3p_dist];
        gene.finalize_coverage(&gene.diff_3p, &mut coverage_3p);

        let mut coverage_pct = vec![0u32; 100];
        gene.finalize_coverage(&gene.diff_percentile, &mut coverage_pct);

        let total_pct_count: u32 = coverage_pct.iter().sum();

        // Skip transcripts with negligible coverage to prevent 1-read amplification spikes
        if total_pct_count < 100 {
            continue;
        }

        active_genes += 1;

        // For normal mode, normalization is per-base coverage.
        // For ends mode, total_pct_count is the number of reads.
        let mean_cov: f64 = if is_ends_mode {
            // In ends mode, we want the 3' profile to represent "probability of having an end at this distance"
            // but scaled so it is comparable across samples.
            (total_pct_count as f64 / 1000.0).max(0.001)
        } else {
            (total_pct_count as f64 / tx_len as f64).max(0.001)
        };

        for (i, &count) in coverage_3p.iter().enumerate() {
            if i >= max_3p_dist {
                break;
            }
            if i as u64 >= tx_len {
                break;
            }

            dist_3p_sums[i] += count as f64 / mean_cov;
            dist_3p_support[i] += 1;
        }

        let norm_pct: f64 = if is_ends_mode {
            // Scale by total reads / 100 so average height is ~1.0
            (total_pct_count as f64 / 100.0).max(0.001)
        } else {
            (total_pct_count as f64 / 100.0).max(0.001)
        };
        for (i, &count) in coverage_pct.iter().enumerate() {
            percentile_sums[i] += count as f64 / norm_pct;
            percentile_support[i] += 1;
        }
    }

    let last_valid = dist_3p_support
        .iter()
        .rposition(|&s| s >= min_support)
        .unwrap_or(0);
    let dynamic_min_support = min_support.max(active_genes / 20); // 5% floor

    let mut dist_3p_means: Vec<f64> = dist_3p_sums
        .into_iter()
        .enumerate()
        .map(|(i, sum)| {
            if dist_3p_support[i] >= dynamic_min_support {
                sum / dist_3p_support[i] as f64
            } else {
                0.0
            }
        })
        .collect();

    // SMOOTHING: Rolling average (window size 50)
    let window = 50;
    if dist_3p_means.len() > window {
        let mut smoothed = Vec::with_capacity(dist_3p_means.len());
        for i in 0..dist_3p_means.len() {
            let start = i.saturating_sub(window / 2);
            let end = (i + window / 2).min(dist_3p_means.len());
            let sum: f64 = dist_3p_means[start..end].iter().sum();
            smoothed.push(sum / (end - start) as f64);
        }
        dist_3p_means = smoothed;
    }

    dist_3p_means.truncate(last_valid + 1);

    AggregatedStats {
        dist_3p_means,
        dist_3p_support,
        percentile_means: percentile_sums
            .into_iter()
            .enumerate()
            .map(|(i, s)| s / percentile_support[i] as f64)
            .collect(),
        percentile_support,
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

    let mut profile_sums = vec![0.0f64; 100];
    let mut total_weight = 0.0f64;
    let mut total_genes = 0usize;

    for gene in index.genes.iter() {
        if gene.total_len < 100 {
            continue;
        }

        let mut coverage_pct = vec![0u32; 100];
        gene.finalize_coverage(&gene.diff_percentile, &mut coverage_pct);

        let total_count: u32 = coverage_pct.iter().sum();
        if total_count < min_support {
            continue;
        }

        // Per-gene normalise: mean bin depth becomes 1.0 (shape only).
        // This is the key step that prevents high-coverage or short transcripts
        // from dominating the aggregate — matches RSeQC's approach.
        let mean_depth = (total_count as f64) / 100.0;

        let weight = if coverage_weighted {
            // Optional: reward well-covered genes, cap to prevent outlier dominance.
            (total_count as f64).min(weight_cap)
        } else {
            // Default (RSeQC-compatible): every gene contributes equally.
            1.0
        };

        for (i, &count) in coverage_pct.iter().enumerate() {
            profile_sums[i] += weight * (count as f64 / mean_depth);
        }
        total_weight += weight;
        total_genes += 1;
    }

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
    let mut accum: Vec<(Vec<f64>, f64, usize)> = BINS
        .iter()
        .map(|_| (vec![0.0f64; 100], 0.0f64, 0usize))
        .collect();

    for gene in index.genes.iter() {
        let tx_len = gene.total_len;

        let Some(bin_idx) = BINS
            .iter()
            .position(|&(_, lo, hi)| tx_len >= lo && tx_len < hi)
        else {
            continue;
        };

        let mut coverage_pct = vec![0u32; 100];
        gene.finalize_coverage(&gene.diff_percentile, &mut coverage_pct);

        let total_count: u32 = coverage_pct.iter().sum();
        if total_count < min_support {
            continue;
        }

        let mean_depth = (total_count as f64) / 100.0;
        let weight = if coverage_weighted {
            (total_count as f64).min(weight_cap)
        } else {
            1.0
        };

        let (sums, tw, n) = &mut accum[bin_idx];
        for (i, &count) in coverage_pct.iter().enumerate() {
            sums[i] += weight * (count as f64 / mean_depth);
        }
        *tw += weight;
        *n += 1;
    }

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
}
