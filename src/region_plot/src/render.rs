use crate::layout::place_reads;
use crate::model::{GeneModel, ReadSegment, RegionPlot, SamplePlotData};
use crate::style::PlotStyle;
use anyhow::{anyhow, Context, Result};
use image::codecs::png::PngEncoder;
use image::{ColorType, ImageEncoder};
use plotters::coord::Shift;
use plotters::prelude::*;
use plotters::style::{register_font, FontStyle};
use std::path::Path;

const EMBEDDED_FONT: &[u8] = include_bytes!("../assets/DejaVuSans.ttf");

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    Png,
    Svg,
}

#[derive(Clone, Debug)]
pub struct PlotOptions {
    pub width: u32,
    pub min_height: u32,
    pub lane_height: u32,
    pub coverage_height: u32,
    pub reference_height: u32,
    pub gene_height: u32,
    pub margin_left: u32,
    pub margin_right: u32,
    pub margin_top: u32,
    pub margin_bottom: u32,
    pub font_family: String,
    pub style: PlotStyle,
}

impl Default for PlotOptions {
    fn default() -> Self {
        Self {
            width: 1400,
            min_height: 500,
            lane_height: 14,
            coverage_height: 80,
            reference_height: 28,
            gene_height: 54,
            margin_left: 96,
            margin_right: 28,
            margin_top: 34,
            margin_bottom: 34,
            font_family: "sans-serif".to_string(),
            style: PlotStyle::default(),
        }
    }
}

pub fn render_png(plot: &RegionPlot, opts: &PlotOptions) -> Result<Vec<u8>> {
    validate_plot(plot)?;
    register_embedded_font(&opts.font_family)?;
    let height = resolved_height(plot, opts);
    let mut rgb = vec![255u8; opts.width as usize * height as usize * 3];

    {
        let root = BitMapBackend::with_buffer(&mut rgb, (opts.width, height)).into_drawing_area();
        draw_region_plot(&root, plot, opts)?;
        root.present().map_err(|e| anyhow!("{e:?}"))?;
    }

    let mut png = Vec::new();
    PngEncoder::new(&mut png)
        .write_image(&rgb, opts.width, height, ColorType::Rgb8)
        .context("failed to encode PNG")?;
    Ok(png)
}

pub fn render_svg(plot: &RegionPlot, opts: &PlotOptions) -> Result<String> {
    validate_plot(plot)?;
    register_embedded_font(&opts.font_family)?;
    let height = resolved_height(plot, opts);
    let mut svg = String::new();

    {
        let root = SVGBackend::with_string(&mut svg, (opts.width, height)).into_drawing_area();
        draw_region_plot(&root, plot, opts)?;
        root.present().map_err(|e| anyhow!("{e:?}"))?;
    }

    Ok(svg)
}

fn register_embedded_font(family: &str) -> Result<()> {
    register_font(family, FontStyle::Normal, EMBEDDED_FONT)
        .map_err(|_| anyhow!("failed to register embedded font"))
}

pub fn render_to_path(plot: &RegionPlot, opts: &PlotOptions, path: impl AsRef<Path>) -> Result<()> {
    let path = path.as_ref();
    let format = match path.extension().and_then(|ext| ext.to_str()) {
        Some(ext) if ext.eq_ignore_ascii_case("png") => OutputFormat::Png,
        Some(ext) if ext.eq_ignore_ascii_case("svg") => OutputFormat::Svg,
        other => return Err(anyhow!("unsupported output extension: {:?}", other)),
    };

    match format {
        OutputFormat::Png => std::fs::write(path, render_png(plot, opts)?)
            .with_context(|| format!("failed to write {}", path.display())),
        OutputFormat::Svg => std::fs::write(path, render_svg(plot, opts)?)
            .with_context(|| format!("failed to write {}", path.display())),
    }
}

fn validate_plot(plot: &RegionPlot) -> Result<()> {
    if plot.end <= plot.start {
        return Err(anyhow!(
            "invalid region {}:{}-{}",
            plot.chrom,
            plot.start,
            plot.end
        ));
    }
    Ok(())
}

fn resolved_height(plot: &RegionPlot, opts: &PlotOptions) -> u32 {
    let read_height: u32 = plot
        .samples
        .iter()
        .map(|sample| {
            let (_, lanes) = place_reads(&sample.reads);
            32 + opts.coverage_height + lanes.max(1) as u32 * opts.lane_height
        })
        .sum();
    let reference = if plot.reference.is_some() {
        opts.reference_height
    } else {
        0
    };
    let genes = if plot.genes.is_empty() {
        0
    } else {
        opts.gene_height
    };

    (opts.margin_top + 34 + reference + genes + read_height + opts.margin_bottom)
        .max(opts.min_height)
}

fn draw_region_plot<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    plot: &RegionPlot,
    opts: &PlotOptions,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.fill(&opts.style.background)?;
    let width = root.dim_in_pixel().0;
    let plot_left = opts.margin_left as i32;
    let plot_right = width as i32 - opts.margin_right as i32;
    let plot_width = (plot_right - plot_left).max(1);

    let mut y = opts.margin_top as i32;
    draw_title(root, plot, opts, y)?;
    y += 34;

    if plot.reference.is_some() {
        draw_reference(root, plot, opts, plot_left, y, plot_width)?;
        y += opts.reference_height as i32;
    }

    if !plot.genes.is_empty() {
        draw_genes(root, &plot.genes, plot, opts, plot_left, y, plot_width)?;
        y += opts.gene_height as i32;
    }

    for sample in &plot.samples {
        y = draw_sample(root, sample, plot, opts, plot_left, y, plot_width)?;
    }

    draw_axis(root, plot, opts, plot_left, y + 8, plot_width)?;
    Ok(())
}

fn draw_title<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    plot: &RegionPlot,
    opts: &PlotOptions,
    y: i32,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let title = format!("{}:{}-{}", plot.chrom, plot.start, plot.end);
    root.draw(&Text::new(
        title,
        (opts.margin_left as i32, y + 18),
        (opts.font_family.as_str(), 22)
            .into_font()
            .color(&opts.style.text),
    ))?;
    Ok(())
}

fn draw_reference<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    plot: &RegionPlot,
    opts: &PlotOptions,
    x0: i32,
    y: i32,
    width: i32,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.draw(&Text::new(
        "Reference",
        (12, y + 18),
        (opts.font_family.as_str(), 14)
            .into_font()
            .color(&opts.style.text),
    ))?;

    let Some(reference) = plot.reference.as_ref() else {
        return Ok(());
    };
    let base_count = reference.len().max(1);
    let max_labels = (width / 14).max(1) as usize;
    let step = (base_count / max_labels).max(1);

    for (idx, base) in reference.iter().enumerate().step_by(step) {
        let x = x0 + ((idx as f64 / base_count as f64) * width as f64) as i32;
        root.draw(&Text::new(
            (*base as char).to_string(),
            (x, y + 18),
            (opts.font_family.as_str(), 12)
                .into_font()
                .color(&opts.style.axis),
        ))?;
    }
    Ok(())
}

fn draw_genes<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    genes: &[GeneModel],
    plot: &RegionPlot,
    opts: &PlotOptions,
    x0: i32,
    y: i32,
    width: i32,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.draw(&Text::new(
        "Genes",
        (12, y + 22),
        (opts.font_family.as_str(), 14)
            .into_font()
            .color(&opts.style.text),
    ))?;

    for (idx, gene) in genes.iter().enumerate() {
        let lane_y = y + 16 + (idx.min(1) as i32 * 18);
        let gx1 = map_pos(gene.start, plot, x0, width);
        let gx2 = map_pos(gene.end, plot, x0, width).max(gx1 + 1);
        root.draw(&PathElement::new(
            vec![(gx1, lane_y), (gx2, lane_y)],
            ShapeStyle::from(&opts.style.gene).stroke_width(2),
        ))?;
        for &(start, end) in &gene.exons {
            let ex1 = map_pos(start, plot, x0, width);
            let ex2 = map_pos(end, plot, x0, width).max(ex1 + 2);
            root.draw(&Rectangle::new(
                [(ex1, lane_y - 5), (ex2, lane_y + 5)],
                opts.style.exon.filled(),
            ))?;
        }
        let label = match gene.strand {
            Some(strand) => format!("{} ({})", gene.name, strand),
            None => gene.name.clone(),
        };
        root.draw(&Text::new(
            label,
            (gx1, lane_y - 8),
            (opts.font_family.as_str(), 11)
                .into_font()
                .color(&opts.style.text),
        ))?;
    }
    Ok(())
}

fn draw_sample<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    sample: &SamplePlotData,
    plot: &RegionPlot,
    opts: &PlotOptions,
    x0: i32,
    y: i32,
    width: i32,
) -> Result<i32>
where
    DB::ErrorType: 'static,
{
    root.draw(&Text::new(
        sample.name.clone(),
        (12, y + 18),
        (opts.font_family.as_str(), 14)
            .into_font()
            .color(&opts.style.text),
    ))?;

    let coverage_top = y + 26;
    draw_coverage(root, sample, plot, opts, x0, coverage_top, width)?;

    let reads_top = coverage_top + opts.coverage_height as i32 + 8;
    let (placed_reads, lanes) = place_reads(&sample.reads);
    for placed in placed_reads {
        let lane_y = reads_top + placed.lane as i32 * opts.lane_height as i32;
        draw_read(root, placed.read, plot, opts, x0, lane_y, width)?;
    }

    Ok(reads_top + lanes.max(1) as i32 * opts.lane_height as i32 + 18)
}

fn draw_coverage<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    sample: &SamplePlotData,
    plot: &RegionPlot,
    opts: &PlotOptions,
    x0: i32,
    y: i32,
    width: i32,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let h = opts.coverage_height as i32;
    root.draw(&PathElement::new(
        vec![(x0, y + h), (x0 + width, y + h)],
        ShapeStyle::from(&opts.style.axis).stroke_width(1),
    ))?;

    let max_depth = sample
        .coverage
        .iter()
        .map(|point| point.depth)
        .max()
        .unwrap_or(1)
        .max(1) as f64;

    for point in &sample.coverage {
        if point.pos < plot.start || point.pos > plot.end {
            continue;
        }
        let x = map_pos(point.pos, plot, x0, width);
        let bar_h = ((point.depth as f64 / max_depth) * h as f64) as i32;
        root.draw(&PathElement::new(
            vec![(x, y + h), (x, y + h - bar_h)],
            ShapeStyle::from(&opts.style.coverage).stroke_width(1),
        ))?;
    }
    Ok(())
}

fn draw_read<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    read: &crate::model::ReadModel,
    plot: &RegionPlot,
    opts: &PlotOptions,
    x0: i32,
    y: i32,
    width: i32,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let color = match read.haplotype {
        Some(1) => opts.style.haplotype_1,
        Some(2) => opts.style.haplotype_2,
        _ if read.is_reverse => opts.style.read_reverse,
        _ => opts.style.read_forward,
    };

    for segment in &read.segments {
        match segment {
            ReadSegment::Match { ref_start, len, .. } => {
                let x1 = map_pos(*ref_start, plot, x0, width);
                let x2 = map_pos(*ref_start + *len, plot, x0, width).max(x1 + 1);
                root.draw(&Rectangle::new([(x1, y), (x2, y + 8)], color.filled()))?;
            }
            ReadSegment::Del { ref_start, len } => {
                let x1 = map_pos(*ref_start, plot, x0, width);
                let x2 = map_pos(*ref_start + *len, plot, x0, width).max(x1 + 1);
                root.draw(&PathElement::new(
                    vec![(x1, y + 4), (x2, y + 4)],
                    ShapeStyle::from(&opts.style.deletion).stroke_width(2),
                ))?;
            }
            ReadSegment::Ins { ref_pos, .. } => {
                let x = map_pos(*ref_pos, plot, x0, width);
                root.draw(&PathElement::new(
                    vec![(x, y - 2), (x, y + 10)],
                    ShapeStyle::from(&opts.style.insertion).stroke_width(2),
                ))?;
            }
            ReadSegment::SoftClip { ref_pos, bases } => {
                let x = map_pos(*ref_pos, plot, x0, width);
                let len = bases.len().max(1) as i32;
                root.draw(&Rectangle::new(
                    [(x, y + 2), (x + len.min(12), y + 6)],
                    ShapeStyle::from(&color.mix(0.45)).filled(),
                ))?;
            }
        }
    }

    for modification in &read.modifications {
        let x = map_pos(modification.ref_pos, plot, x0, width);
        root.draw(&Circle::new(
            (x, y - 1),
            2,
            opts.style.modification.filled(),
        ))?;
    }

    Ok(())
}

fn draw_axis<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    plot: &RegionPlot,
    opts: &PlotOptions,
    x0: i32,
    y: i32,
    width: i32,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.draw(&PathElement::new(
        vec![(x0, y), (x0 + width, y)],
        ShapeStyle::from(&opts.style.axis).stroke_width(1),
    ))?;
    let ticks = 5;
    for idx in 0..=ticks {
        let frac = idx as f64 / ticks as f64;
        let pos = plot.start + (plot.span() as f64 * frac).round() as i64;
        let x = x0 + (width as f64 * frac) as i32;
        root.draw(&PathElement::new(
            vec![(x, y), (x, y + 5)],
            ShapeStyle::from(&opts.style.axis).stroke_width(1),
        ))?;
        root.draw(&Text::new(
            pos.to_string(),
            (x - 18, y + 20),
            (opts.font_family.as_str(), 11)
                .into_font()
                .color(&opts.style.axis),
        ))?;
    }
    Ok(())
}

fn map_pos(pos: i64, plot: &RegionPlot, x0: i32, width: i32) -> i32 {
    let frac = ((pos - plot.start) as f64 / plot.span() as f64).clamp(0.0, 1.0);
    x0 + (frac * width as f64) as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{CoveragePoint, ReadModel};

    fn sample_plot() -> RegionPlot {
        RegionPlot {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            reference: Some(b"ACGTACGTACGT".to_vec()),
            genes: vec![GeneModel {
                name: "GENE1".to_string(),
                start: 110,
                end: 190,
                strand: Some('+'),
                exons: vec![(120, 140), (160, 180)],
            }],
            samples: vec![SamplePlotData {
                name: "sample".to_string(),
                coverage: vec![
                    CoveragePoint { pos: 100, depth: 3 },
                    CoveragePoint { pos: 150, depth: 9 },
                ],
                reads: vec![ReadModel {
                    name: "read".to_string(),
                    start: 110,
                    end: 180,
                    is_reverse: false,
                    mapq: 60,
                    segments: vec![ReadSegment::Match {
                        ref_start: 110,
                        len: 70,
                        query_start: 0,
                    }],
                    bases: None,
                    qualities: None,
                    haplotype: Some(1),
                    modifications: Vec::new(),
                }],
            }],
        }
    }

    #[test]
    fn renders_png_bytes() {
        let bytes = render_png(&sample_plot(), &PlotOptions::default()).unwrap();
        assert!(bytes.starts_with(b"\x89PNG"));
    }

    #[test]
    fn renders_svg_text() {
        let svg = render_svg(&sample_plot(), &PlotOptions::default()).unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("GENE1"));
    }
}
