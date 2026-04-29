use anyhow::Result;
use kuva::prelude::*;
use std::collections::HashMap;

/// Generates a combined Multi-Sample plot for Gene Body Coverage (RSeQC-classic).
pub fn generate_gene_body_plot(
    all_data: &HashMap<String, Vec<f64>>,
    output_path: &str,
) -> Result<()> {
    let mut plots: Vec<Plot> = Vec::new();

    // Sort keys for consistent legend order and color assignment
    let mut samples: Vec<_> = all_data.keys().collect();
    samples.sort();

    for sample_id in samples {
        let y_values = all_data.get(sample_id).unwrap();
        let data: Vec<(f64, f64)> = y_values
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as f64 + 1.0, v))
            .collect();

        let line = LinePlot::new()
            .with_data(data)
            .with_legend(sample_id.clone())
            .with_line_style(LineStyle::Solid);

        plots.push(line.into());
    }

    let layout = Layout::auto_from_plots(&plots)
        .with_title("Gene Body Coverage (RSeQC-Classic)")
        .with_x_label("Percentile (5' -> 3')")
        .with_y_label("Percentage Coverage (0-100%)");

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
    output_path: &str,
) -> Result<()> {
    let mut plots: Vec<Plot> = Vec::new();

    // Sort samples for consistent legend order and colors
    let mut samples: Vec<_> = all_data.keys().collect();
    samples.sort();

    for sample_id in samples {
        let y_values = all_data.get(sample_id).unwrap();
        let data: Vec<(f64, f64)> = y_values
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as f64, v))
            .collect();

        let line = LinePlot::new()
            .with_data(data)
            .with_legend(sample_id.clone())
            .with_line_style(LineStyle::Solid);

        plots.push(line.into());
    }

    let layout = Layout::auto_from_plots(&plots)
        .with_title("3' Distance Bias Profile (Combined)")
        .with_x_label("Distance from 3' end (bp)")
        .with_y_label("Summed Normalized Coverage");

    let svg = render_to_svg(plots, layout);
    std::fs::write(output_path, svg)?;

    Ok(())
}

/// Generates a single profile plot for 3' Distance Bias.
pub fn generate_3p_dist_plot(sample_id: &str, data: &[f64], output_path: &str) -> Result<()> {
    let mut all_data = HashMap::new();
    all_data.insert(sample_id.to_string(), data.to_vec());
    generate_multi_3p_dist_plot(&all_data, output_path)
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
