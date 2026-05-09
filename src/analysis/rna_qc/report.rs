use super::config::RnaQcConfig;
use crate::analysis::qc::RnaSeqQcSummary;
use crate::analysis::snapshot::{resolve_snapshot_region_with_flank, run_snapshot, SnapshotConfig};
use crate::io::annotation::AnnotationFormat;
use crate::models::Gene;
use crate::stats::plotting::{
    generate_gene_body_plot, generate_multi_3p_dist_plot, generate_multi_inner_distance_plot,
    generate_multi_raw_3p_plot, generate_rna_qc_summary_svg, generate_stratified_gene_body_plot,
    render_rna_qc_terminal_summary, PlotMetadata,
};
use crate::stats::{
    write_3p_wide_tsv, write_classic_wide_format, write_feature_counts_tsv,
    write_gene_profiles_tsv, AggregatedStats,
};
use anyhow::{Context, Result};
use serde::Serialize;
use std::collections::HashMap;
use std::fs;
use std::io::Write;

#[derive(Debug, Serialize)]
struct RnaSampleSummary {
    sample: String,
    total_reads: u64,
    active_genes: usize,
    active_3p_genes: usize,
    aligned_qc_reads: u64,
    mtdna_reads: u64,
    mtdna_fraction: f64,
    rdna_reads: u64,
    rdna_fraction: f64,
    requested_pairs: usize,
    informative_pairs: usize,
    inferred_strandness: String,
    fr_count: usize,
    rf_count: usize,
    other_orientation_count: usize,
    inner_distance_mean: Option<f64>,
    inner_distance_median: Option<f64>,
    exonic_reads: u64,
    intronic_reads: u64,
    flanking_reads: u64,
    intergenic_reads: u64,
}

pub(crate) fn write_sample_reports(
    config: &RnaQcConfig,
    sample_name: &str,
    stats: &AggregatedStats,
) -> Result<()> {
    if config.write_gene_profiles
        && config
            .analysis
            .contains(&crate::analysis::types::AnalysisType::ThreePrime)
    {
        let tsv_start = std::time::Instant::now();
        let profiles_path = format!("{}.{}.gene_profiles.tsv", config.output, sample_name);
        write_gene_profiles_tsv(&profiles_path, &stats.gene_qc, config.compact_gene_profiles)?;
        println!(
            "Gene profiles written to: {} (took: {:?})",
            profiles_path,
            tsv_start.elapsed()
        );
    }

    if config.write_counts {
        let counts_path = format!("{}.{}.gene_counts.tsv", config.output, sample_name);
        write_feature_counts_tsv(&counts_path, &stats.gene_qc)?;
        println!("Feature-count-like matrix written to: {}", counts_path);
    }

    Ok(())
}

pub(crate) fn write_sample_visual_reports(
    config: &RnaQcConfig,
    bam_path: &str,
    sample_name: &str,
    stats: &AggregatedStats,
    qc: &RnaSeqQcSummary,
    read_dist_counts: &[u64; 12],
    genes: &[Gene],
    annotation_format: AnnotationFormat,
) -> Result<()> {
    let term_start = std::time::Instant::now();
    print_terminal_summary(sample_name, stats, qc, read_dist_counts)?;
    println!("  - Terminal summary rendering took: {:?}", term_start.elapsed());
    write_sample_summary_json(config, sample_name, stats, qc, read_dist_counts)?;

    if !config.no_plot {
        let svg_path = format!("{}.{}.rna.qc_summary.svg", config.output, sample_name);
        generate_rna_qc_summary_svg(sample_name, stats, qc, read_dist_counts, &svg_path)?;
        println!("  - RNA QC summary figure written to: {}", svg_path);
    }

    if config.snap_qc || !config.snap_genes.is_empty() {
        write_snapshots(
            config,
            bam_path,
            sample_name,
            genes,
            annotation_format,
            stats,
        )?;
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
    )?;

    if config.stratify_length && !stats.stratified_percentile_normalized.is_empty() {
        generate_stratified_gene_body_plot(
            &stats.stratified_percentile_normalized,
            &format!(
                "{}.{}.geneBodyCoverage.stratified.svg",
                config.output, sample_name
            ),
        )?;
    }

    Ok(())
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

fn print_terminal_summary(
    sample_name: &str,
    stats: &AggregatedStats,
    qc: &RnaSeqQcSummary,
    read_dist_counts: &[u64; 12],
) -> Result<()> {
    let rendered = render_rna_qc_terminal_summary(sample_name, stats, qc, read_dist_counts)?;
    println!("\n{}", rendered);
    println!("RNA QC summary");
    println!("--------------");
    println!("Reads passing QC:        {}", qc.aligned_qc_reads);
    println!(
        "mtDNA fraction:          {:.1}%",
        qc.mtdna_fraction() * 100.0
    );
    println!(
        "Stranded reverse frac:   {:.1}%",
        stranded_reverse_percent(qc)
    );
    let read_dist = read_distribution_percentages(read_dist_counts);
    println!(
        "Read distribution:       exon {:.1}% ({}) | intron {:.1}% ({}) | flank {:.1}% ({}) | intergenic {:.1}% ({})",
        read_dist.exonic, read_dist.exonic_count,
        read_dist.intronic, read_dist.intronic_count,
        read_dist.flanking, read_dist.flanking_count,
        read_dist.intergenic, read_dist.intergenic_count
    );
    println!("Strandness:              {}", qc.inferred_strandness());
    println!(
        "Inner distance median:   {}",
        qc.inner_distance_median()
            .map(|v| format!("{v:.0} bp"))
            .unwrap_or_else(|| "NA".to_string())
    );
    println!("Gene-body genes:         {}", stats.active_genes);
    println!("3' profile genes:        {}", stats.active_3p_genes);
    if stats.unknown_chrom_reads > 0 {
        println!("Unknown chrom reads:     {}", stats.unknown_chrom_reads);
    }
    println!("Total tags assigned:     {}", stats.total_tags);
    Ok(())
}

fn stranded_reverse_percent(qc: &RnaSeqQcSummary) -> f64 {
    let stranded_total = qc.stranded_forward_count + qc.stranded_reverse_count;
    if stranded_total == 0 {
        0.0
    } else {
        percent_from_counts(qc.stranded_reverse_count as u64, stranded_total as f64)
    }
}

struct ReadDistributionPercentages {
    exonic: f64,
    intronic: f64,
    flanking: f64,
    intergenic: f64,
    exonic_count: u64,
    intronic_count: u64,
    flanking_count: u64,
    intergenic_count: u64,
}

fn read_distribution_percentages(read_dist_counts: &[u64; 12]) -> ReadDistributionPercentages {
    let exonic = read_dist_counts[crate::analysis::read_distribution::RegionType::CdsExon as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::Utr5Exon as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::Utr3Exon as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::Exon as usize];
    let intronic =
        read_dist_counts[crate::analysis::read_distribution::RegionType::Intron as usize];
    let flanking = read_dist_counts
        [crate::analysis::read_distribution::RegionType::TssUp1kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TssUp5kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TssUp10kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TesDown1kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TesDown5kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TesDown10kb as usize];
    let intergenic =
        read_dist_counts[crate::analysis::read_distribution::RegionType::Intergenic as usize];
    let total = (exonic + intronic + flanking + intergenic).max(1) as f64;
    ReadDistributionPercentages {
        exonic: percent_from_counts(exonic, total),
        intronic: percent_from_counts(intronic, total),
        flanking: percent_from_counts(flanking, total),
        intergenic: percent_from_counts(intergenic, total),
        exonic_count: exonic,
        intronic_count: intronic,
        flanking_count: flanking,
        intergenic_count: intergenic,
    }
}

fn percent_from_counts(count: u64, total: f64) -> f64 {
    if total <= 0.0 {
        0.0
    } else {
        (count as f64 * 100.0 / total).clamp(0.0, 100.0)
    }
}

fn write_sample_summary_json(
    config: &RnaQcConfig,
    sample_name: &str,
    stats: &AggregatedStats,
    qc: &RnaSeqQcSummary,
    read_dist_counts: &[u64; 12],
) -> Result<()> {
    let summary = build_sample_summary(sample_name, stats, qc, read_dist_counts);
    crate::analysis::report::write_summary_json(
        &format!("{}.{}.rna.summary.json", config.output, sample_name),
        "rna",
        sample_name,
        &summary,
    )?;
    println!(
        "  - RNA summary JSON written to: {}.{}.rna.summary.json",
        config.output, sample_name
    );
    Ok(())
}

fn build_sample_summary(
    sample_name: &str,
    stats: &AggregatedStats,
    qc: &RnaSeqQcSummary,
    read_dist_counts: &[u64; 12],
) -> RnaSampleSummary {
    let exonic = read_dist_counts[crate::analysis::read_distribution::RegionType::CdsExon as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::Utr5Exon as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::Utr3Exon as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::Exon as usize];
    let intronic =
        read_dist_counts[crate::analysis::read_distribution::RegionType::Intron as usize];
    let flanking = read_dist_counts
        [crate::analysis::read_distribution::RegionType::TssUp1kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TssUp5kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TssUp10kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TesDown1kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TesDown5kb as usize]
        + read_dist_counts[crate::analysis::read_distribution::RegionType::TesDown10kb as usize];
    let intergenic =
        read_dist_counts[crate::analysis::read_distribution::RegionType::Intergenic as usize];

    RnaSampleSummary {
        sample: sample_name.to_string(),
        total_reads: stats.total_reads,
        active_genes: stats.active_genes,
        active_3p_genes: stats.active_3p_genes,
        aligned_qc_reads: qc.aligned_qc_reads,
        mtdna_reads: qc.mtdna_reads,
        mtdna_fraction: qc.mtdna_fraction(),
        rdna_reads: qc.rdna_reads,
        rdna_fraction: qc.rdna_fraction(),
        requested_pairs: qc.requested_pairs,
        informative_pairs: qc.informative_pairs,
        inferred_strandness: qc.inferred_strandness().to_string(),
        fr_count: qc.fr_count,
        rf_count: qc.rf_count,
        other_orientation_count: qc.other_orientation_count,
        inner_distance_mean: qc.inner_distance_mean(),
        inner_distance_median: qc.inner_distance_median(),
        exonic_reads: exonic,
        intronic_reads: intronic,
        flanking_reads: flanking,
        intergenic_reads: intergenic,
    }
}

fn write_snapshots(
    config: &RnaQcConfig,
    bam_path: &str,
    sample_name: &str,
    genes: &[Gene],
    annotation_format: AnnotationFormat,
    stats: &AggregatedStats,
) -> Result<()> {
    let snap_dir = format!("{}.{}.rna_snapshots", config.output, sample_name);
    fs::create_dir_all(&snap_dir)
        .with_context(|| format!("could not create snapshot directory {snap_dir}"))?;
    let manifest_path = format!("{}.{}.rna_snapshots.tsv", config.output, sample_name);
    let mut manifest = fs::File::create(&manifest_path)
        .with_context(|| format!("could not create snapshot manifest {manifest_path}"))?;
    writeln!(
        manifest,
        "gene\tchrom\tstart\tend\tstrand\tsnapshot_path\treason"
    )?;

    for target in snapshot_targets(config, genes, stats) {
        match resolve_snapshot_region_with_flank(
            &target.gene_name,
            Some(&config.annotation),
            annotation_format,
            Some(config.snap_flank),
        ) {
            Ok(region) => {
                let safe_name = sanitize_snapshot_label(&target.gene_name);
                let snapshot_path = format!("{snap_dir}/{safe_name}.png");
                let snapshot_cfg = SnapshotConfig {
                    bam_path: bam_path.to_string(),
                    bai_path: None,
                    region: region.clone(),
                    annotation_path: Some(config.annotation.clone()),
                    annotation_format,
                    reference_path: None,
                    output_path: snapshot_path.clone(),
                    mapq_threshold: config.mapq,
                    max_reads: config.snap_max_reads,
                    width: 1400,
                    min_height: 500,
                    show_reference: false,
                    show_genes: true,
                    show_reference_base_track: false,
                    show_sample_base_track: false,
                    squash: false,
                };
                run_snapshot(&snapshot_cfg).with_context(|| {
                    format!("failed to render snapshot for {}", target.gene_name)
                })?;
                writeln!(
                    manifest,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    target.gene_name,
                    region.chrom,
                    region.start + 1,
                    region.end,
                    target.strand,
                    snapshot_path,
                    target.reason
                )?;
            }
            Err(err) => {
                writeln!(manifest, "{}\t\t\t\t\t\t{}", target.gene_name, err)?;
            }
        }
    }

    println!("  - Snapshot manifest written to: {}", manifest_path);
    println!("  - Snapshot directory written to: {}", snap_dir);
    Ok(())
}

struct SnapshotTarget {
    gene_name: String,
    strand: char,
    reason: String,
}

fn snapshot_targets(
    config: &RnaQcConfig,
    genes: &[Gene],
    stats: &AggregatedStats,
) -> Vec<SnapshotTarget> {
    let mut targets = Vec::new();
    for name in &config.snap_genes {
        if let Some(gene) = find_gene_by_name(genes, name) {
            targets.push(SnapshotTarget {
                gene_name: gene.name.clone().unwrap_or_else(|| name.clone()),
                strand: gene.representative.strand,
                reason: "custom_gene".to_string(),
            });
        } else {
            targets.push(SnapshotTarget {
                gene_name: name.clone(),
                strand: '+',
                reason: "custom_gene_missing".to_string(),
            });
        }
    }

    if config.snap_qc {
        for name in default_snapshot_genes() {
            if targets
                .iter()
                .any(|t| t.gene_name.eq_ignore_ascii_case(&name))
            {
                continue;
            }
            if let Some(gene) = find_gene_by_name(genes, &name) {
                targets.push(SnapshotTarget {
                    gene_name: gene.name.clone().unwrap_or(name.clone()),
                    strand: gene.representative.strand,
                    reason: "default_marker".to_string(),
                });
            }
        }
        if stats.active_genes > 0 {
            if let Some(long_gene) = genes
                .iter()
                .filter(|g| g.name.as_deref().is_some())
                .max_by_key(|g| g.representative.total_length)
            {
                let long_name = long_gene.name.clone().unwrap();
                if !targets
                    .iter()
                    .any(|t| t.gene_name.eq_ignore_ascii_case(&long_name))
                {
                    targets.push(SnapshotTarget {
                        gene_name: long_name,
                        strand: long_gene.representative.strand,
                        reason: "long_gene".to_string(),
                    });
                }
            }
        }
    }

    targets
}

fn default_snapshot_genes() -> Vec<String> {
    vec![
        "GAPDH".to_string(),
        "ACTB".to_string(),
        "MALAT1".to_string(),
        "RPLP0".to_string(),
        "RPS18".to_string(),
        "MT-CO1".to_string(),
        "MT-ND1".to_string(),
    ]
}

fn find_gene_by_name<'a>(genes: &'a [Gene], name: &str) -> Option<&'a Gene> {
    let needle = name.trim().to_lowercase();
    genes.iter().find(|gene| {
        gene.name
            .as_deref()
            .map(|n| n.to_lowercase() == needle)
            .unwrap_or_else(|| gene.id.to_lowercase() == needle)
    })
}

fn sanitize_snapshot_label(label: &str) -> String {
    let mut out = String::new();
    for ch in label.chars() {
        if ch.is_ascii_alphanumeric() || ch == '_' || ch == '-' {
            out.push(ch);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        "snapshot".to_string()
    } else {
        out
    }
}
