use super::config::RnaQcConfig;
use crate::analysis::index::AnnotationIndex;
use crate::stats::{
    aggregate_genes, aggregate_rseqc_classic, aggregate_rseqc_stratified, AggregatedStats,
    ThreePrimeParams,
};

pub(crate) fn aggregate_sample(
    index: &AnnotationIndex,
    config: &RnaQcConfig,
    state: &crate::analysis::rna_qc::state::RnaWorkerState,
    ref_juncs: &crate::analysis::splice_junction::ReferenceJunctions,
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

    let mut stats = aggregate_genes(
        index,
        config.min_support,
        &three_prime_params,
        state,
        false,
    );

    let (gene_body_raw, gene_body_normalized) = aggregate_rseqc_classic(
        index,
        false,
        config.min_support as u32,
        config.coverage_weighted,
    );
    stats.percentile_means = gene_body_raw;
    stats.percentile_normalized = gene_body_normalized;

    if config.stratify_length {
        stats.stratified_percentile_normalized =
            aggregate_rseqc_stratified(index, config.min_support as u32, config.coverage_weighted);
    }

    stats.splice_junctions = Some(crate::analysis::splice_junction::compute_metrics(
        &state.splice_junctions,
        ref_juncs,
    ));

    stats
}
