use crate::model::RegionPlot;
use crate::scene::{build_scene, Scene, VisualElement};
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
    pub ruler_height: u32,
    pub reference_height: u32,
    pub base_track_height: u32,
    pub gene_height: u32,
    pub sample_label_height: u32,
    pub read_track_gap: u32,
    pub sample_gap: u32,
    pub margin_left: u32,
    pub margin_right: u32,
    pub margin_top: u32,
    pub margin_bottom: u32,
    pub min_alt_af: f32,
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
            ruler_height: 42,
            reference_height: 28,
            base_track_height: 28,
            gene_height: 54,
            sample_label_height: 26,
            read_track_gap: 8,
            sample_gap: 18,
            margin_left: 96,
            margin_right: 28,
            margin_top: 34,
            margin_bottom: 34,
            min_alt_af: 0.20,
            font_family: "sans-serif".to_string(),
            style: PlotStyle::default(),
        }
    }
}

pub fn render_png(plot: &RegionPlot, opts: &PlotOptions) -> Result<Vec<u8>> {
    validate_plot(plot)?;
    register_embedded_font(&opts.font_family)?;
    let scene = build_scene(plot, opts);
    let mut rgb = vec![255u8; scene.width as usize * scene.height as usize * 3];

    {
        let root =
            BitMapBackend::with_buffer(&mut rgb, (scene.width, scene.height)).into_drawing_area();
        draw_scene(&root, &scene, opts)?;
        root.present().map_err(|e| anyhow!("{e:?}"))?;
    }

    let mut png = Vec::new();
    PngEncoder::new(&mut png)
        .write_image(&rgb, scene.width, scene.height, ColorType::Rgb8)
        .context("failed to encode PNG")?;
    Ok(png)
}

pub fn render_svg(plot: &RegionPlot, opts: &PlotOptions) -> Result<String> {
    validate_plot(plot)?;
    register_embedded_font(&opts.font_family)?;
    let scene = build_scene(plot, opts);
    let mut svg = String::new();

    {
        let root =
            SVGBackend::with_string(&mut svg, (scene.width, scene.height)).into_drawing_area();
        draw_scene(&root, &scene, opts)?;
        root.present().map_err(|e| anyhow!("{e:?}"))?;
    }

    Ok(svg)
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

fn register_embedded_font(family: &str) -> Result<()> {
    register_font(family, FontStyle::Normal, EMBEDDED_FONT)
        .map_err(|_| anyhow!("failed to register embedded font"))
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

fn draw_scene<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    scene: &Scene,
    opts: &PlotOptions,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    for element in &scene.elements {
        match element {
            VisualElement::Rect {
                x1,
                y1,
                x2,
                y2,
                color,
            } => {
                let color: RGBColor = (*color).into();
                root.draw(&Rectangle::new(
                    [((*x1) as i32, (*y1) as i32), ((*x2) as i32, (*y2) as i32)],
                    color.filled(),
                ))?;
            }
            VisualElement::Line {
                x1,
                y1,
                x2,
                y2,
                color,
                width,
            } => {
                let color: RGBColor = (*color).into();
                root.draw(&PathElement::new(
                    vec![((*x1) as i32, (*y1) as i32), ((*x2) as i32, (*y2) as i32)],
                    ShapeStyle::from(&color).stroke_width(*width),
                ))?;
            }
            VisualElement::Text {
                x,
                y,
                text,
                color,
                size,
            } => {
                let color: RGBColor = (*color).into();
                root.draw(&Text::new(
                    text.clone(),
                    ((*x) as i32, (*y) as i32),
                    (opts.font_family.as_str(), *size).into_font().color(&color),
                ))?;
            }
            VisualElement::Circle {
                x,
                y,
                radius,
                color,
            } => {
                let color: RGBColor = (*color).into();
                root.draw(&Circle::new(
                    ((*x) as i32, (*y) as i32),
                    *radius as i32,
                    color.filled(),
                ))?;
            }
            VisualElement::Triangle { x, y, size, color } => {
                let color: RGBColor = (*color).into();
                let s = *size;
                root.draw(&Polygon::new(
                    vec![
                        ((*x) as i32, (*y - s) as i32),
                        ((*x - s) as i32, (*y + s) as i32),
                        ((*x + s) as i32, (*y + s) as i32),
                    ],
                    color.filled(),
                ))?;
            }
            VisualElement::Polygon { points, color } => {
                let color: RGBColor = (*color).into();
                root.draw(&Polygon::new(
                    points
                        .iter()
                        .map(|(x, y)| ((*x) as i32, (*y) as i32))
                        .collect::<Vec<_>>(),
                    color.filled(),
                ))?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{
        BasePileup, CoveragePoint, GeneModel, ReadModel, ReadSegment, SamplePlotData,
    };

    fn sample_plot() -> RegionPlot {
        RegionPlot {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            reference: Some(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec()),
            genes: vec![GeneModel {
                name: "GENE1".to_string(),
                start: 110,
                end: 190,
                strand: Some('+'),
                exons: vec![(120, 140), (160, 180)],
            }],
            samples: vec![SamplePlotData {
                name: "sample".to_string(),
                pileup: (0..100)
                    .map(|idx| BasePileup {
                        a: if idx % 4 == 0 { 8 } else { 0 },
                        c: if idx % 4 == 1 { 8 } else { 0 },
                        g: if idx % 4 == 2 { 8 } else { 0 },
                        t: if idx % 4 == 3 { 8 } else { 0 },
                        total: 8 + (idx % 9) as u32,
                        ..BasePileup::default()
                    })
                    .collect(),
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
                    segments: vec![
                        ReadSegment::Match {
                            ref_start: 110,
                            len: 25,
                            query_start: 0,
                        },
                        ReadSegment::Mismatch {
                            ref_start: 135,
                            len: 1,
                            query_start: 25,
                            base: b'T',
                        },
                        ReadSegment::Ins {
                            ref_pos: 150,
                            bases: b"AC".to_vec(),
                        },
                        ReadSegment::Del {
                            ref_start: 160,
                            len: 5,
                        },
                    ],
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

    #[test]
    fn scene_has_rich_tracks() {
        let scene = build_scene(&sample_plot(), &PlotOptions::default());
        assert!(scene.elements.len() > 100);
    }
}
