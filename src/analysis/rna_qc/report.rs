use super::config::RnaQcConfig;
use crate::stats::plotting::{
    generate_gene_body_plot, generate_multi_3p_dist_plot, generate_multi_inner_distance_plot,
    generate_multi_raw_3p_plot, PlotMetadata,
};
use crate::stats::{
    write_3p_wide_tsv, write_classic_wide_format, write_feature_counts_tsv,
    write_gene_profiles_tsv, AggregatedStats,
};
use anyhow::Result;
use std::collections::HashMap;

pub(crate) fn write_sample_reports(
    config: &RnaQcConfig,
    sample_name: &str,
    stats: &AggregatedStats,
) -> Result<()> {
    if config
        .analysis
        .contains(&crate::analysis::types::AnalysisType::ThreePrime)
    {
        let profiles_path = format!("{}.{}.gene_profiles.tsv", config.output, sample_name);
        write_gene_profiles_tsv(&profiles_path, &stats.gene_qc)?;
        println!("Gene profiles written to: {}", profiles_path);
    }

    if config.write_counts {
        let counts_path = format!("{}.{}.gene_counts.tsv", config.output, sample_name);
        write_feature_counts_tsv(&counts_path, &stats.gene_qc)?;
        println!("Feature-count-like matrix written to: {}", counts_path);
    }

    Ok(())
}

pub(crate) fn write_sample_gene_body_plot(
    config: &RnaQcConfig,
    sample_name: &str,
    stats: &AggregatedStats,
) -> Result<()> {
    if config.no_plot
        || !config
            .analysis
            .contains(&crate::analysis::types::AnalysisType::GeneBody)
    {
        return Ok(());
    }

    let data = HashMap::from([(sample_name.to_string(), stats.percentile_means.clone())]);
    let metadata = HashMap::from([(
        sample_name.to_string(),
        PlotMetadata {
            total_reads: stats.total_reads,
            active_genes: stats.active_genes,
            active_3p_genes: stats.active_3p_genes,
        },
    )]);

    generate_gene_body_plot(
        &data,
        &metadata,
        &format!("{}.{}.geneBodyCoverage.svg", config.output, sample_name),
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn write_multi_sample_outputs(
    config: &RnaQcConfig,
    classic_all: &HashMap<String, Vec<f64>>,
    classic_percent_all: &HashMap<String, Vec<f64>>,
    dist_3p_all: &HashMap<String, Vec<f64>>,
    dist_3p_support_all: &HashMap<String, Vec<usize>>,
    dist_3p_raw_all: &HashMap<String, Vec<f64>>,
    dist_3p_raw_analyzed_all: &HashMap<String, Vec<f64>>,
    inner_distance_all: &HashMap<String, Vec<f64>>,
    plot_metadata_all: &HashMap<String, PlotMetadata>,
) -> Result<()> {
    if config
        .analysis
        .contains(&crate::analysis::types::AnalysisType::GeneBody)
        && !classic_all.is_empty()
    {
        let classic_path = format!("{}.geneBodyCoverage.txt", config.output);
        write_classic_wide_format(&classic_path, classic_all)?;
        println!(
            "\nMultiQC-compatible raw aggregate written to: {}",
            classic_path
        );

        let classic_pct_path = format!("{}.geneBodyCoverage.percent.txt", config.output);
        write_classic_wide_format(&classic_pct_path, classic_percent_all)?;
        println!(
            "Normalized percentage aggregate written to: {}",
            classic_pct_path
        );

        if !config.no_plot {
            generate_gene_body_plot(
                classic_all,
                plot_metadata_all,
                &format!("{}.geneBodyCoverage.svg", config.output),
            )?;
            println!(
                "Gene body coverage plot written to: {}.geneBodyCoverage.svg",
                config.output
            );
        }
    }

    if config
        .analysis
        .contains(&crate::analysis::types::AnalysisType::ThreePrime)
        && !dist_3p_all.is_empty()
    {
        let dist_3p_path = format!("{}.3p_anchor_normalized.tsv", config.output);
        write_3p_wide_tsv(
            &dist_3p_path,
            dist_3p_all,
            config.three_prime_bin_size,
            config.max_3p_dist,
        )?;
        println!("3' anchor-normalized profile written to: {}", dist_3p_path);

        let dist_3p_support_path = format!("{}.3p_support.tsv", config.output);
        write_3p_wide_tsv(
            &dist_3p_support_path,
            dist_3p_support_all,
            config.three_prime_bin_size,
            config.max_3p_dist,
        )?;
        println!("3' support profile written to: {}", dist_3p_support_path);

        let dist_3p_raw_analyzed_path = format!("{}.3p_raw_analyzed.tsv", config.output);
        write_3p_wide_tsv(
            &dist_3p_raw_analyzed_path,
            dist_3p_raw_analyzed_all,
            config.three_prime_bin_size,
            config.max_3p_dist,
        )?;
        println!(
            "3' raw analyzed profile written to: {}",
            dist_3p_raw_analyzed_path
        );

        let dist_3p_raw_all_path = format!("{}.3p_raw_all.tsv", config.output);
        write_3p_wide_tsv(
            &dist_3p_raw_all_path,
            dist_3p_raw_all,
            config.three_prime_bin_size,
            config.max_3p_dist,
        )?;
        println!("3' raw all profile written to: {}", dist_3p_raw_all_path);

        if !config.no_plot {
            generate_multi_3p_dist_plot(
                dist_3p_all,
                plot_metadata_all,
                &format!("{}.3p_dist.svg", config.output),
                config.three_prime_bin_size,
            )?;
            generate_multi_raw_3p_plot(
                dist_3p_raw_analyzed_all,
                plot_metadata_all,
                &format!("{}.3p_raw.svg", config.output),
                config.three_prime_bin_size,
            )?;
        }
    }

    if config
        .analysis
        .contains(&crate::analysis::types::AnalysisType::Qc)
        && !inner_distance_all.is_empty()
        && !config.no_plot
    {
        generate_multi_inner_distance_plot(
            inner_distance_all,
            &format!("{}.inner_distance.svg", config.output),
            -200,
        )?;
        println!(
            "Inner distance profile written to: {}.inner_distance.svg",
            config.output
        );
    }

    Ok(())
}
