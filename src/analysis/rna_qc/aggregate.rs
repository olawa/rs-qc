use super::config::RnaQcConfig;
use crate::analysis::index::AnnotationIndex;
use crate::stats::{aggregate_genes, AggregatedStats, ThreePrimeParams};

pub(crate) fn aggregate_sample(
    index: &AnnotationIndex,
    config: &RnaQcConfig,
    total_reads: u64,
) -> AggregatedStats {
    let three_prime_params = ThreePrimeParams {
        normalization_bp: config.normalization_bp,
        max_3p_dist: config.max_3p_dist,
        bin_size: config.three_prime_bin_size,
        min_anchor_count: config.three_prime_min_anchor_count,
        min_anchor_mean: config.three_prime_min_anchor_mean,
        min_anchor_nonzero_bins: config.three_prime_min_anchor_nonzero_bins,
        max_ratio: config.three_prime_max_ratio,
    };

    aggregate_genes(
        index,
        config.min_support,
        &three_prime_params,
        total_reads,
        false,
    )
}
