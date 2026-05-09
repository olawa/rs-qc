use crate::analysis::qc::RnaSeqQcSummary;
use anyhow::Result;
use kuva::backend::svg::SvgBackend;
use kuva::prelude::*;
use std::collections::HashMap;

pub struct PlotMetadata {
    pub total_reads: u64,
    pub active_genes: usize,
    pub active_3p_genes: usize,
}

/// Generates a combined Multi-Sample plot for Gene Body Coverage (RSeQC-classic).
pub fn generate_gene_body_plot(
    all_data: &HashMap<String, Vec<f64>>,
    metadata: &HashMap<String, PlotMetadata>,
    output_path: &str,
) -> Result<()> {
    let mut plots: Vec<Plot> = Vec::new();
    let mut subtitle_parts = Vec::new();

    // Sort keys for consistent legend order and color assignment
    let mut samples: Vec<_> = all_data.keys().collect();
    samples.sort();

    for sample_id in samples {
        let y_values = all_data.get(sample_id).unwrap();
        let max_val = y_values.iter().fold(0.0f64, |a, &b| a.max(b)).max(0.001);
        let data: Vec<(f64, f64)> = y_values
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as f64 + 1.0, (v / max_val) * 100.0))
            .collect();

        let line = LinePlot::new()
            .with_data(data)
            .with_legend(sample_id.clone())
            .with_line_style(LineStyle::Solid);

        plots.push(line.into());

        if let Some(m) = metadata.get(sample_id) {
            subtitle_parts.push(format!(
                "{}: {} reads, {} genes",
                sample_id, m.total_reads, m.active_genes
            ));
        }
    }

    let title = if subtitle_parts.is_empty() {
        "Gene Body Coverage (RSeQC-Classic)".to_string()
    } else {
        format!(
            "Gene Body Coverage (RSeQC-Classic) [{}]",
            subtitle_parts.join(" | ")
        )
    };

    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_x_label("Percentile (5' -> 3')")
        .with_y_label("Percentage Coverage (0-100% of max)");

    let svg = render_to_svg(plots, layout);
    std::fs::write(output_path, svg)?;

    Ok(())
}

/// Generates the length-stratified gene body coverage plot with a fixed, visually
/// distinct colour per length class. The colour order (blue → teal → orange → red)
/// encodes increasing transcript length, making the progression of 5' bias
/// immediately obvious.
///
/// Recognised keys (as emitted by `aggregate_rseqc_stratified`):
///
/// | Key                  | Colour        |
/// |----------------------|---------------|
/// | `short (<1.5kb)`     | `#2196F3` (blue)  |
/// | `medium (1.5-5kb)`   | `#009688` (teal)  |
/// | `long (5-10kb)`      | `#FF9800` (orange)|
/// | `very long (>10kb)`  | `#F44336` (red)   |
///
/// Any key not in the table falls back to grey (`#9E9E9E`).
pub fn generate_stratified_gene_body_plot(
    strat_data: &HashMap<String, Vec<f64>>,
    output_path: &str,
) -> Result<()> {
    // Fixed colour per length class — must match the labels in aggregate_rseqc_stratified.
    let colour_map: &[(&str, &str)] = &[
        ("short (<1.5kb)", "#2196F3"),    // blue
        ("medium (1.5-5kb)", "#009688"),  // teal
        ("long (5-10kb)", "#FF9800"),     // orange
        ("very long (>10kb)", "#F44336"), // red
    ];
    let fallback_colour = "#9E9E9E"; // grey for unexpected keys

    // Preserve biological order rather than alphabetical sort.
    let ordered_keys: Vec<&str> = colour_map
        .iter()
        .map(|&(k, _)| k)
        .filter(|k| strat_data.contains_key(*k))
        .collect();

    let mut plots: Vec<Plot> = Vec::new();

    for &key in &ordered_keys {
        let colour = colour_map
            .iter()
            .find(|&&(k, _)| k == key)
            .map(|&(_, c)| c)
            .unwrap_or(fallback_colour);

        let y_values = strat_data.get(key).unwrap();
        let data: Vec<(f64, f64)> = y_values
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as f64 + 1.0, v))
            .collect();

        let line = LinePlot::new()
            .with_data(data)
            .with_legend(key.to_string())
            .with_color(colour)
            .with_line_style(LineStyle::Solid);

        plots.push(line.into());
    }

    let layout = Layout::auto_from_plots(&plots)
        .with_title("Gene Body Coverage by Transcript Length")
        .with_x_label("Percentile (5' -> 3')")
        .with_y_label("Percentage Coverage (0-100%)");

    let svg = render_to_svg(plots, layout);
    std::fs::write(output_path, svg)?;

    Ok(())
}

/// Generates a combined Multi-Sample plot for 3' Distance Bias.
pub fn generate_multi_3p_dist_plot(
    all_data: &HashMap<String, Vec<f64>>,
    metadata: &HashMap<String, PlotMetadata>,
    output_path: &str,
    bin_size: usize,
) -> Result<()> {
    let mut plots: Vec<Plot> = Vec::new();
    let mut subtitle_parts = Vec::new();

    // Sort samples for consistent legend order and colors
    let mut samples: Vec<_> = all_data.keys().collect();
    samples.sort();

    for sample_id in samples {
        let y_values = all_data.get(sample_id).unwrap();
        let data: Vec<(f64, f64)> = y_values
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as f64 * bin_size as f64, v))
            .collect();

        let line = LinePlot::new()
            .with_data(data)
            .with_legend(sample_id.clone())
            .with_line_style(LineStyle::Solid);

        plots.push(line.into());

        if let Some(m) = metadata.get(sample_id) {
            let active_label = if m.active_3p_genes > 0 {
                m.active_3p_genes
            } else {
                m.active_genes
            };
            subtitle_parts.push(format!(
                "{}: {} reads, {} genes",
                sample_id, m.total_reads, active_label
            ));
        }
    }

    // Add vertical markers at 5kb and 10kb
    for &x in &[5000.0, 10000.0] {
        plots.push(
            LinePlot::new()
                .with_data(vec![(x, 0.0), (x, 2.0)])
                .with_color("#888")
                .with_line_style(LineStyle::Dashed)
                .with_stroke_width(1.0)
                .into(),
        );
    }

    let title = if subtitle_parts.is_empty() {
        "3' Distance Bias Profile (Combined)".to_string()
    } else {
        format!(
            "3' Distance Bias Profile (Combined) [{}]",
            subtitle_parts.join(" | ")
        )
    };

    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_x_label("Distance from 3' end (bp)")
        .with_y_label("Mean normalized coverage relative to the 3' anchor window");

    let svg = render_to_svg(plots, layout);
    std::fs::write(output_path, svg)?;

    Ok(())
}

/// Generates a single profile plot for 3' Distance Bias.
#[allow(dead_code)]
pub fn generate_3p_dist_plot(
    sample_id: &str,
    data: &[f64],
    metadata: &PlotMetadata,
    output_path: &str,
    bin_size: usize,
) -> Result<()> {
    let mut all_data = HashMap::new();
    all_data.insert(sample_id.to_string(), data.to_vec());
    let mut all_meta = HashMap::new();
    all_meta.insert(
        sample_id.to_string(),
        PlotMetadata {
            total_reads: metadata.total_reads,
            active_genes: metadata.active_genes,
            active_3p_genes: metadata.active_3p_genes,
        },
    );
    generate_multi_3p_dist_plot(&all_data, &all_meta, output_path, bin_size)
}

/// Generates a combined Multi-Sample plot for raw 3' distance distribution.
pub fn generate_multi_raw_3p_plot(
    all_data: &HashMap<String, Vec<f64>>,
    metadata: &HashMap<String, PlotMetadata>,
    output_path: &str,
    bin_size: usize,
) -> Result<()> {
    let mut plots: Vec<Plot> = Vec::new();
    let mut subtitle_parts = Vec::new();

    let mut samples: Vec<_> = all_data.keys().collect();
    samples.sort();

    for sample_id in samples {
        let y_values = all_data.get(sample_id).unwrap();
        let data: Vec<(f64, f64)> = y_values
            .iter()
            .enumerate()
            .map(|(i, &v)| ((i * bin_size) as f64, v))
            .collect();

        let line = LinePlot::new()
            .with_data(data)
            .with_legend(sample_id.clone())
            .with_line_style(LineStyle::Solid);

        plots.push(line.into());

        if let Some(m) = metadata.get(sample_id) {
            subtitle_parts.push(format!(
                "{}: {} reads, {} genes",
                sample_id, m.total_reads, m.active_genes
            ));
        }
    }

    let title = if subtitle_parts.is_empty() {
        "Raw 3' Distance Profile (Expression Size)".to_string()
    } else {
        format!(
            "Raw 3' Distance Profile (Expression Size) [{}]",
            subtitle_parts.join(" | ")
        )
    };

    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_x_label("Distance from 3' end (bp)")
        .with_y_label("Total raw counts (bases)");

    let svg = render_to_svg(plots, layout);
    std::fs::write(output_path, svg)?;

    Ok(())
}

/// Generates a combined Multi-Sample plot for inner distance around the fragment center.
pub fn generate_multi_inner_distance_plot(
    all_data: &HashMap<String, Vec<f64>>,
    output_path: &str,
    min_dist: i64,
) -> Result<()> {
    let mut plots: Vec<Plot> = Vec::new();

    let mut samples: Vec<_> = all_data.keys().collect();
    samples.sort();

    for sample_id in samples {
        let y_values = all_data.get(sample_id).unwrap();
        let data: Vec<(f64, f64)> = y_values
            .iter()
            .enumerate()
            .map(|(i, &v)| (min_dist as f64 + i as f64, v))
            .collect();

        let line = LinePlot::new()
            .with_data(data)
            .with_legend(sample_id.clone())
            .with_line_style(LineStyle::Solid);

        plots.push(line.into());
    }

    let layout = Layout::auto_from_plots(&plots)
        .with_title("Inner Distance Distribution")
        .with_x_label("Inner distance (bp)")
        .with_y_label("Sampled pair count");

    let svg = render_to_svg(plots, layout);
    std::fs::write(output_path, svg)?;

    Ok(())
}

pub fn generate_rna_qc_summary_svg(
    sample_name: &str,
    stats: &crate::stats::AggregatedStats,
    qc: &RnaSeqQcSummary,
    read_dist_counts: &[u64; 12],
    output_path: &str,
) -> Result<()> {
    let figure = build_rna_qc_summary_figure(sample_name, stats, qc, read_dist_counts, true, false)
        .with_figure_size(1600.0, 980.0)
        .with_title(format!("RNA QC summary for {}", sample_name))
        .with_title_size(24)
        .with_labels();
    let scene = figure.render();
    std::fs::write(output_path, SvgBackend.render_scene(&scene))?;
    Ok(())
}

pub fn render_rna_qc_terminal_summary(
    sample_name: &str,
    stats: &crate::stats::AggregatedStats,
    qc: &RnaSeqQcSummary,
    read_dist_counts: &[u64; 12],
) -> Result<String> {
    let figure = build_rna_qc_summary_figure(sample_name, stats, qc, read_dist_counts, false, true)
        .with_figure_size(1180.0, 740.0)
        .with_title(format!("RNA QC summary for {}", sample_name))
        .with_title_size(14);
    let scene = figure.render();
    let cols = std::env::var("COLUMNS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map(|n| n.clamp(60, 80))
        .unwrap_or(80);
    let rows = std::env::var("LINES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map(|n| n.clamp(32, 48))
        .unwrap_or(40);
    Ok(TerminalBackend::new(cols, rows).render_scene(&scene))
}

fn build_rna_qc_summary_figure(
    sample_name: &str,
    stats: &crate::stats::AggregatedStats,
    qc: &RnaSeqQcSummary,
    read_dist_counts: &[u64; 12],
    include_raw_profile: bool,
    compact_labels: bool,
) -> Figure {
    let gene_body_max = stats
        .percentile_means
        .iter()
        .fold(0.0_f64, |acc, &v| acc.max(v))
        .max(0.001);
    let gene_body_plot = LinePlot::new()
        .with_data(
            stats
                .percentile_means
                .iter()
                .enumerate()
                .map(|(i, &v)| (i as f64 + 1.0, (v / gene_body_max) * 100.0)),
        )
        .with_color("steelblue")
        .with_stroke_width(2.0);

    let three_prime_plot = LinePlot::new()
        .with_data(
            stats
                .dist_3p_means
                .iter()
                .enumerate()
                .map(|(i, &v)| (i as f64 * stats.bin_size as f64, v)),
        )
        .with_color("darkorange")
        .with_stroke_width(2.0);

    let mut three_prime_plots = vec![Plot::Line(three_prime_plot)];
    for &x in &[5000.0, 10000.0] {
        three_prime_plots.push(
            LinePlot::new()
                .with_data(vec![(x, 0.0), (x, 1.5)])
                .with_color("#aaa")
                .with_line_style(LineStyle::Dashed)
                .with_stroke_width(1.0)
                .into(),
        );
    }

    let raw_three_prime_plot = LinePlot::new()
        .with_data(
            stats
                .dist_3p_sums_raw_analyzed
                .iter()
                .enumerate()
                .map(|(i, &v)| (i as f64 * stats.bin_size as f64, v)),
        )
        .with_color("seagreen")
        .with_stroke_width(2.0);

    let inner_distance_plot = LinePlot::new()
        .with_data(
            qc.clipped_inner_distance_series(-200, 200)
                .into_iter()
                .enumerate()
                .map(|(i, v)| (-200.0 + i as f64, v)),
        )
        .with_color("crimson")
        .with_stroke_width(2.0);

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
    let read_total = (exonic + intronic + flanking + intergenic).max(1) as f64;
    let read_distribution_plot = category_bar_plot(vec![
        ("Exonic", percent_from_counts(exonic, read_total), "#4e79a7"),
        (
            "Intronic",
            percent_from_counts(intronic, read_total),
            "#59a14f",
        ),
        (
            "Flanking",
            percent_from_counts(flanking, read_total),
            "#f28e2b",
        ),
        (
            "Intergenic",
            percent_from_counts(intergenic, read_total),
            "#b07aa1",
        ),
        (
            "mtDNA",
            percent_from_fraction(qc.mtdna_fraction()),
            "#e15759",
        ),
        ("Rev strand", stranded_reverse_percent(qc), "#9c755f"),
    ]);

    let qc_plot = category_bar_plot(vec![
        ("mtDNA", qc.mtdna_fraction() * 100.0, "#e15759"),
        ("rDNA", qc.rdna_fraction() * 100.0, "#76b7b2"),
        ("FR", qc.fr_fraction() * 100.0, "#59a14f"),
        ("RF", qc.rf_fraction() * 100.0, "#4e79a7"),
        ("Other", qc.other_fraction() * 100.0, "#edc948"),
    ]);

    let mut plots: Vec<Vec<Plot>> = vec![
        vec![Plot::Line(gene_body_plot)],
        three_prime_plots,
        vec![Plot::Line(inner_distance_plot)],
        vec![Plot::Bar(read_distribution_plot)],
    ];
    let mut layouts = vec![
        Layout::auto_from_plots(&plots[0])
            .with_title(if compact_labels {
                "Gene body"
            } else {
                "Gene body"
            })
            .with_x_label(if compact_labels {
                ""
            } else {
                "5' -> 3' percentile"
            })
            .with_y_label(if compact_labels {
                ""
            } else {
                "Relative coverage"
            }),
        Layout::auto_from_plots(&plots[1])
            .with_title(if compact_labels {
                "3' profile"
            } else {
                "3' anchor profile"
            })
            .with_x_label(if compact_labels {
                ""
            } else {
                "Distance from 3' end (bp)"
            })
            .with_y_label(if compact_labels {
                ""
            } else {
                "Normalized coverage"
            }),
        Layout::auto_from_plots(&plots[2])
            .with_title(if compact_labels {
                "Inner dist"
            } else {
                "Inner distance"
            })
            .with_x_label(if compact_labels { "" } else { "bp" })
            .with_y_label(if compact_labels { "" } else { "Count" }),
        Layout::auto_from_plots(&plots[3])
            .with_title(if compact_labels {
                "Read dist"
            } else {
                "Read distribution"
            })
            .with_x_label(if compact_labels {
                ""
            } else {
                "Region / metric"
            })
            .with_y_label(if compact_labels { "" } else { "Share (%)" }),
    ];

    let read_layout = std::mem::replace(&mut layouts[3], Layout::auto_from_plots(&plots[3]))
        .with_y_axis_min(0.0)
        .with_y_axis_max(100.0)
        .with_y_tick_format(TickFormat::Fixed(0));
    layouts[3] = read_layout;

    if include_raw_profile {
        plots.insert(2, vec![Plot::Line(raw_three_prime_plot)]);
        layouts.insert(
            2,
            Layout::auto_from_plots(&plots[2])
                .with_title(if compact_labels {
                    "Raw 3'"
                } else {
                    "Raw 3' profile"
                })
                .with_x_label(if compact_labels {
                    ""
                } else {
                    "Distance from 3' end (bp)"
                })
                .with_y_label(if compact_labels { "" } else { "Raw counts" }),
        );
    }

    if include_raw_profile {
        plots.push(vec![Plot::Bar(qc_plot)]);
        layouts.push(
            Layout::auto_from_plots(plots.last().unwrap())
                .with_title("QC fractions")
                .with_x_label("")
                .with_y_label(""),
        );
    }

    let rows = if include_raw_profile { 2 } else { 2 };
    let cols = if include_raw_profile { 3 } else { 2 };
    let mut figure = Figure::new(rows, cols)
        .with_plots(plots)
        .with_layouts(layouts);
    if !compact_labels {
        figure = figure.with_labels();
    }
    figure = figure.with_title(format!("RNA QC summary for {}", sample_name));
    figure
}

fn category_bar_plot(categories: Vec<(&str, f64, &str)>) -> BarPlot {
    let mut plot = BarPlot::new();
    for (label, value, color) in categories {
        plot = plot.with_group(label, vec![(value, color)]);
    }
    plot
}

fn percent_from_counts(count: u64, total: f64) -> f64 {
    if total <= 0.0 {
        0.0
    } else {
        (count as f64 * 100.0 / total).clamp(0.0, 100.0)
    }
}

fn percent_from_fraction(fraction: f64) -> f64 {
    (fraction * 100.0).clamp(0.0, 100.0)
}

fn stranded_reverse_percent(qc: &RnaSeqQcSummary) -> f64 {
    let stranded_total = qc.stranded_forward_count + qc.stranded_reverse_count;
    if stranded_total == 0 {
        0.0
    } else {
        percent_from_counts(qc.stranded_reverse_count as u64, stranded_total as f64)
    }
}
