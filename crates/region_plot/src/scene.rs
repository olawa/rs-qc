use crate::layout::{place_genes, place_reads};
use crate::model::{
    BasePileup, GeneModel, MarkerType, ReadModel, ReadSegment, RegionPlot, SamplePlotData,
};
use crate::render::PlotOptions;
use plotters::style::RGBColor;
use std::collections::HashSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Rgb(pub u8, pub u8, pub u8);

impl From<RGBColor> for Rgb {
    fn from(value: RGBColor) -> Self {
        Self(value.0, value.1, value.2)
    }
}

impl From<Rgb> for RGBColor {
    fn from(value: Rgb) -> Self {
        RGBColor(value.0, value.1, value.2)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum VisualElement {
    Rect {
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        color: Rgb,
    },
    Line {
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        color: Rgb,
        width: u32,
    },
    Text {
        x: f64,
        y: f64,
        text: String,
        color: Rgb,
        size: u32,
    },
    Circle {
        x: f64,
        y: f64,
        radius: f64,
        color: Rgb,
    },
    Triangle {
        x: f64,
        y: f64,
        size: f64,
        color: Rgb,
    },
    Polygon {
        points: Vec<(f64, f64)>,
        color: Rgb,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct Scene {
    pub width: u32,
    pub height: u32,
    pub elements: Vec<VisualElement>,
}

#[derive(Clone, Copy)]
struct PlotGeom {
    x0: f64,
    width: f64,
}

impl PlotGeom {
    fn x(&self, pos: i64, plot: &RegionPlot) -> f64 {
        let frac = ((pos - plot.start) as f64 / plot.span() as f64).clamp(0.0, 1.0);
        self.x0 + frac * self.width
    }

    fn base_x2(&self, pos: i64, plot: &RegionPlot) -> f64 {
        self.x(pos + 1, plot).max(self.x(pos, plot) + 0.75)
    }
}

fn draw_background_markers(
    elements: &mut Vec<VisualElement>,
    plot: &RegionPlot,
    geom: PlotGeom,
    height: f64,
) {
    for marker in &plot.markers {
        if matches!(marker.marker_type, MarkerType::RegionOfInterest) {
            let start = marker.pos;
            let end = marker.end_pos.unwrap_or(start + 1);
            if start >= plot.end || end <= plot.start {
                continue;
            }
            let x1 = geom.x(start, plot);
            let x2 = geom.x(end, plot).max(x1 + 1.0);
            elements.push(VisualElement::Rect {
                x1,
                y1: 0.0,
                x2,
                y2: height,
                color: Rgb(254, 243, 199), // Soft amber background highlight
            });
        }
    }
}

fn draw_foreground_markers(
    elements: &mut Vec<VisualElement>,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    height: f64,
) {
    let mut label_y_offset = 0.0;
    for marker in &plot.markers {
        let x = geom.x(marker.pos, plot);
        match marker.marker_type {
            MarkerType::Variant => {
                if marker.pos >= plot.start && marker.pos < plot.end {
                    elements.push(VisualElement::Line {
                        x1: x,
                        y1: opts.margin_top as f64,
                        x2: x,
                        y2: height - opts.margin_bottom as f64,
                        color: Rgb(220, 38, 38), // Red
                        width: 2,
                    });

                    let pin_y = opts.margin_top as f64 + 4.0;
                    elements.push(VisualElement::Circle {
                        x,
                        y: pin_y,
                        radius: 4.5,
                        color: Rgb(220, 38, 38),
                    });

                    elements.push(VisualElement::Text {
                        x: (x + 6.0).min(opts.width as f64 - 150.0),
                        y: pin_y + 4.0 + label_y_offset,
                        text: marker.label.clone(),
                        color: Rgb(220, 38, 38),
                        size: 11,
                    });
                    label_y_offset = (label_y_offset + 12.0) % 36.0;
                }
            }
            MarkerType::StructuralVariant => {
                let start = marker.pos;
                let end = marker.end_pos.unwrap_or(start + 1);
                if start >= plot.end || end <= plot.start {
                    continue;
                }
                let x1 = geom.x(start, plot);
                let x2 = geom.x(end, plot).max(x1 + 1.0);
                let bracket_y = opts.margin_top as f64 + 20.0;

                elements.push(VisualElement::Line {
                    x1,
                    y1: bracket_y,
                    x2,
                    y2: bracket_y,
                    color: Rgb(37, 99, 235), // Blue
                    width: 3,
                });
                elements.push(VisualElement::Line {
                    x1,
                    y1: bracket_y - 4.0,
                    x2: x1,
                    y2: bracket_y + 4.0,
                    color: Rgb(37, 99, 235),
                    width: 2,
                });
                elements.push(VisualElement::Line {
                    x1: x2,
                    y1: bracket_y - 4.0,
                    x2: x2,
                    y2: bracket_y + 4.0,
                    color: Rgb(37, 99, 235),
                    width: 2,
                });

                for bx in &[x1, x2] {
                    elements.push(VisualElement::Line {
                        x1: *bx,
                        y1: bracket_y + 4.0,
                        x2: *bx,
                        y2: height - opts.margin_bottom as f64,
                        color: Rgb(147, 197, 253), // Soft light blue boundary
                        width: 1,
                    });
                }

                let text_x = x1 + (x2 - x1) * 0.5 - (marker.label.len() as f64 * 3.0);
                elements.push(VisualElement::Text {
                    x: text_x.clamp(geom.x0, opts.width as f64 - 150.0),
                    y: bracket_y - 6.0,
                    text: marker.label.clone(),
                    color: Rgb(37, 99, 235),
                    size: 11,
                });
            }
            MarkerType::RegionOfInterest => {}
        }
    }
}

pub fn build_scene(plot: &RegionPlot, opts: &PlotOptions) -> Scene {
    let height = resolved_height(plot, opts);
    let geom = PlotGeom {
        x0: opts.margin_left as f64,
        width: (opts.width - opts.margin_left - opts.margin_right).max(1) as f64,
    };
    let mut elements = vec![VisualElement::Rect {
        x1: 0.0,
        y1: 0.0,
        x2: opts.width as f64,
        y2: height as f64,
        color: opts.style.background.into(),
    }];

    draw_background_markers(&mut elements, plot, geom, height as f64);

    let mut y = opts.margin_top as f64;
    elements.push(VisualElement::Text {
        x: geom.x0,
        y: y + 20.0,
        text: format!("{}:{}-{}", plot.chrom, plot.start, plot.end),
        color: opts.style.text.into(),
        size: 22,
    });
    y += 34.0;

    draw_ruler(&mut elements, plot, opts, geom, y);
    y += opts.ruler_height as f64;

    if !plot.genes.is_empty() {
        draw_genes(&mut elements, &plot.genes, plot, opts, geom, y);
        y += gene_track_height(plot, opts) as f64;
    }

    if opts.show_reference_base_track {
        if let Some(reference) = &plot.reference {
            draw_reference_strip(&mut elements, reference, plot, opts, geom, y);
            y += opts.reference_height as f64;
        }
    }

    for sample in &plot.samples {
        draw_sample(&mut elements, sample, plot, opts, geom, y);
        y += sample_height(sample, opts) as f64;
    }

    draw_foreground_markers(&mut elements, plot, opts, geom, height as f64);

    Scene {
        width: opts.width,
        height,
        elements,
    }
}

pub(crate) fn resolved_height(plot: &RegionPlot, opts: &PlotOptions) -> u32 {
    let sample_heights: u32 = plot
        .samples
        .iter()
        .map(|sample| sample_height(sample, opts))
        .sum();
    let gene_h = gene_track_height(plot, opts);
    let reference_h = if plot.reference.is_some() && opts.show_reference_base_track {
        opts.reference_height
    } else {
        0
    };
    (opts.margin_top
        + 34
        + opts.ruler_height
        + gene_h
        + reference_h
        + sample_heights
        + opts.margin_bottom)
        .max(opts.min_height)
}

fn sample_height(sample: &SamplePlotData, opts: &PlotOptions) -> u32 {
    let (_, lanes) = place_reads(&sample.reads, opts.squash);
    let lane_h = if opts.squash { 3 } else { 14 };
    let gap = if opts.squash { 0 } else { 2 };
    let base_track_h = if opts.show_sample_base_track {
        opts.reference_height
    } else {
        0
    };
    opts.sample_label_height
        + opts.coverage_height
        + base_track_h
        + if base_track_h > 0 {
            opts.read_track_gap
        } else {
            opts.read_track_gap / 2
        }
        + (lanes as u32 * (lane_h + gap))
        + opts.sample_gap
}

fn gene_track_height(plot: &RegionPlot, _opts: &PlotOptions) -> u32 {
    if plot.genes.is_empty() {
        0
    } else {
        let (_, lanes) = place_genes(&plot.genes);
        (lanes.max(1) as u32 * 16) + 10
    }
}

fn draw_ruler(
    elements: &mut Vec<VisualElement>,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    elements.push(VisualElement::Line {
        x1: geom.x0,
        y1: y + 16.0,
        x2: geom.x0 + geom.width,
        y2: y + 16.0,
        color: opts.style.axis.into(),
        width: 1,
    });
    for idx in 0..=5 {
        let frac = idx as f64 / 5.0;
        let pos = plot.start + (plot.span() as f64 * frac).round() as i64;
        let x = geom.x0 + geom.width * frac;
        elements.push(VisualElement::Line {
            x1: x,
            y1: y + 11.0,
            x2: x,
            y2: y + 21.0,
            color: opts.style.axis.into(),
            width: 1,
        });
        elements.push(VisualElement::Text {
            x: x - 24.0,
            y: y + 34.0,
            text: pos.to_string(),
            color: opts.style.axis.into(),
            size: 11,
        });
    }
}

fn draw_genes(
    elements: &mut Vec<VisualElement>,
    genes: &[GeneModel],
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    elements.push(VisualElement::Text {
        x: 12.0,
        y: y + 16.0,
        text: "Genes".to_string(),
        color: opts.style.text.into(),
        size: 13,
    });
    let (placed_genes, _) = place_genes(genes);
    let mut labeled_names = HashSet::new();
    for placed in placed_genes {
        let gene = placed.gene;
        // Draw intron lines with strand arrows
        let lane_y = y + 14.0 + placed.lane as f64 * 16.0;
        let x1 = geom.x(gene.start, plot);
        let x2 = geom.x(gene.end, plot).max(x1 + 1.0);
        let (gene_color, exon_color) = gene_colors(gene, opts);
        elements.push(VisualElement::Line {
            x1,
            y1: lane_y,
            x2,
            y2: lane_y,
            color: gene_color.into(),
            width: 1,
        });

        // Add strand arrows
        if let Some(strand) = gene.strand {
            let px_per_base = geom.width / plot.span() as f64;
            let arrow_spacing = (100.0 / px_per_base).max(20.0);
            let mut arrow_pos = (gene.start as f64 / arrow_spacing).ceil() * arrow_spacing;
            while arrow_pos < gene.end as f64 {
                if arrow_pos > gene.start as f64 {
                    let ax = geom.x(arrow_pos as i64, plot);
                    if ax > x1 + 5.0 && ax < x2 - 5.0 {
                        let dx1 = if strand == '+' { -3.0 } else { 3.0 };
                        elements.push(VisualElement::Line {
                            x1: ax + dx1,
                            y1: lane_y - 3.0,
                            x2: ax,
                            y2: lane_y,
                            color: gene_color.into(),
                            width: 1,
                        });
                        elements.push(VisualElement::Line {
                            x1: ax + dx1,
                            y1: lane_y + 3.0,
                            x2: ax,
                            y2: lane_y,
                            color: gene_color.into(),
                            width: 1,
                        });
                    }
                }
                arrow_pos += arrow_spacing;
            }
        }

        for &(start, end) in &gene.exons {
            let ex1 = geom.x(start, plot);
            let ex2 = geom.x(end, plot).max(ex1 + 2.0);
            elements.push(VisualElement::Rect {
                x1: ex1,
                y1: lane_y - 4.0,
                x2: ex2,
                y2: lane_y + 4.0,
                color: exon_color.into(),
            });
        }
        if x2 - x1 > 16.0 && labeled_names.insert(gene.name.clone()) {
            let text_x = x1.max(geom.x0 + 4.0);
            elements.push(VisualElement::Text {
                x: text_x,
                y: lane_y - 10.0,
                text: gene.name.clone(),
                color: opts.style.text.into(),
                size: 11,
            });
        }
    }
}

enum BaseSource<'a> {
    Ref(&'a [u8]),
    Pileup(&'a [BasePileup]),
}

impl<'a> BaseSource<'a> {
    fn len(&self) -> usize {
        match self {
            BaseSource::Ref(r) => r.len(),
            BaseSource::Pileup(p) => p.len(),
        }
    }

    fn get_counts(&self, idx: usize) -> [u32; 5] {
        match self {
            BaseSource::Ref(r) => {
                let mut c = [0; 5];
                if let Some(&b) = r.get(idx) {
                    match b.to_ascii_uppercase() {
                        b'A' => c[0] = 1,
                        b'C' => c[1] = 1,
                        b'G' => c[2] = 1,
                        b'T' => c[3] = 1,
                        _ => c[4] = 1,
                    }
                }
                c
            }
            BaseSource::Pileup(p) => {
                let mut c = [0; 5];
                if let Some(pile) = p.get(idx) {
                    c[0] = pile.a;
                    c[1] = pile.c;
                    c[2] = pile.g;
                    c[3] = pile.t;
                    c[4] = pile.n;
                }
                c
            }
        }
    }

    fn get_range_counts(&self, start: usize, end: usize) -> [u32; 5] {
        let mut total = [0; 5];
        for i in start..end {
            let c = self.get_counts(i);
            for j in 0..5 {
                total[j] += c[j];
            }
        }
        total
    }

    fn dominant_base_in_range(&self, start: usize, end: usize) -> u8 {
        let counts = self.get_range_counts(start, end);
        let (idx, _) = counts
            .iter()
            .enumerate()
            .max_by_key(|(_, &count)| count)
            .unwrap_or((4, &0));
        [b'A', b'C', b'G', b'T', b'N'][idx]
    }
}

fn draw_base_track_internal(
    elements: &mut Vec<VisualElement>,
    source: BaseSource,
    label: &str,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
    height: f64,
    show_allele_freq: bool,
    margin: f64,
) {
    if !label.is_empty() {
        elements.push(VisualElement::Text {
            x: 12.0,
            y: y + 18.0,
            text: label.to_string(),
            color: opts.style.text.into(),
            size: 13,
        });
    }

    let px_per_base = geom.width / plot.span() as f64;
    let draw_text = px_per_base > 10.0;

    if px_per_base >= 1.0 {
        for idx in 0..source.len() {
            let pos = plot.start + idx as i64;
            if pos >= plot.end {
                break;
            }

            let counts = source.get_counts(idx);
            let total: u32 = counts.iter().sum::<u32>().max(1);

            let mut lanes = Vec::new();
            for (i, &base) in [b'A', b'C', b'G', b'T'].iter().enumerate() {
                if counts[i] > 0 {
                    lanes.push((base, counts[i] as f32 / total as f32));
                }
            }
            lanes.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

            let x1 = geom.x(pos, plot);
            let x2 = geom.base_x2(pos, plot);
            let mut current_y = y + margin;
            let track_h = height - 2.0 * margin;

            if !show_allele_freq || lanes.len() <= 1 {
                let base = source.dominant_base_in_range(idx, idx + 1);
                if base != b'N' {
                    elements.push(VisualElement::Rect {
                        x1,
                        y1: current_y,
                        x2,
                        y2: current_y + track_h,
                        color: base_color(base).into(),
                    });
                    if draw_text {
                        let text_size = (track_h * 0.85).clamp(8.0, 16.0) as u32;
                        let text_x = x1 + (x2 - x1) * 0.5 - (text_size as f64 * 0.3);
                        let text_y = current_y + track_h * 0.5 + (text_size as f64 * 0.35);
                        elements.push(VisualElement::Text {
                            x: text_x,
                            y: text_y,
                            text: (base as char).to_string(),
                            color: Rgb(255, 255, 255),
                            size: text_size,
                        });
                    }
                }
            } else {
                for (base, frac) in lanes {
                    let h = track_h * frac as f64;
                    let y1 = current_y;
                    let y2 = y1 + h;
                    current_y = y2;
                    elements.push(VisualElement::Rect {
                        x1,
                        y1,
                        x2,
                        y2,
                        color: base_color(base).into(),
                    });
                    if draw_text && h > 8.0 {
                        let text_size = (h * 0.85).clamp(8.0, 16.0) as u32;
                        let text_x = x1 + (x2 - x1) * 0.5 - (text_size as f64 * 0.3);
                        let text_y = y1 + h * 0.5 + (text_size as f64 * 0.35);
                        elements.push(VisualElement::Text {
                            x: text_x,
                            y: text_y,
                            text: (base as char).to_string(),
                            color: Rgb(255, 255, 255),
                            size: text_size,
                        });
                    }
                }
            }
        }
    } else {
        // Aggregated view (Zoom out) - Always use dominant base to avoid "scrambled" look
        let bases_per_px = 1.0 / px_per_base;
        let track_h = height - 2.0 * margin;
        for x_px in 0..geom.width as i32 {
            let b_start = (x_px as f64 * bases_per_px).floor() as usize;
            let b_end = ((x_px + 1) as f64 * bases_per_px).ceil() as usize;
            let b_end = b_end.min(source.len());
            if b_start >= b_end {
                continue;
            }

            let base = source.dominant_base_in_range(b_start, b_end);
            if base == b'N' {
                continue;
            }

            let x1 = geom.x0 + x_px as f64;
            elements.push(VisualElement::Rect {
                x1,
                y1: y + margin,
                x2: x1 + 1.0,
                y2: y + margin + track_h,
                color: base_color(base).into(),
            });
        }
    }
}

fn gene_colors(gene: &GeneModel, opts: &PlotOptions) -> (RGBColor, RGBColor) {
    if gene.strand == Some('-') {
        (opts.style.gene_reverse, opts.style.exon_reverse)
    } else {
        (opts.style.gene, opts.style.exon)
    }
}

fn draw_sample(
    elements: &mut Vec<VisualElement>,
    sample: &SamplePlotData,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    elements.push(VisualElement::Text {
        x: 12.0,
        y: y + 18.0,
        text: sample.name.clone(),
        color: opts.style.text.into(),
        size: 13,
    });
    let coverage_y = y + opts.sample_label_height as f64;
    draw_coverage(elements, sample, plot, opts, geom, coverage_y);
    let base_y = coverage_y + opts.coverage_height as f64;
    let read_y = if opts.show_sample_base_track {
        draw_sample_base_track(elements, sample, plot, opts, geom, base_y);
        base_y + opts.reference_height as f64 + opts.read_track_gap as f64
    } else {
        base_y + (opts.read_track_gap / 2) as f64
    };
    draw_reads(elements, sample, plot, opts, geom, read_y);
}

fn draw_coverage(
    elements: &mut Vec<VisualElement>,
    sample: &SamplePlotData,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    let h = opts.coverage_height as f64;
    elements.push(VisualElement::Line {
        x1: geom.x0,
        y1: y + h,
        x2: geom.x0 + geom.width,
        y2: y + h,
        color: opts.style.coverage.into(),
        width: 1,
    });
    let depths = coverage_depths(sample, plot.start);
    let max_depth = depths
        .iter()
        .map(|(_, depth)| *depth)
        .max()
        .unwrap_or(1)
        .max(1) as f64;
    for (pos, depth) in depths {
        if pos < plot.start || pos >= plot.end {
            continue;
        }
        let x1 = geom.x(pos, plot);
        let bar_h = (depth as f64 / max_depth) * h;
        elements.push(VisualElement::Rect {
            x1,
            y1: y + h - bar_h,
            x2: geom.base_x2(pos, plot).max(x1 + 1.0),
            y2: y + h,
            color: opts.style.coverage.into(),
        });
    }
}

fn draw_sample_base_track(
    elements: &mut Vec<VisualElement>,
    _sample: &SamplePlotData,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    let reference = match plot.reference.as_ref() {
        Some(r) => r,
        None => return,
    };
    draw_base_track_internal(
        elements,
        BaseSource::Ref(reference),
        "",
        plot,
        opts,
        geom,
        y,
        opts.reference_height as f64,
        false,
        4.0,
    );
}

fn draw_reference_strip(
    elements: &mut Vec<VisualElement>,
    reference: &[u8],
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    draw_base_track_internal(
        elements,
        BaseSource::Ref(reference),
        "Reference",
        plot,
        opts,
        geom,
        y,
        opts.reference_height as f64,
        false,
        4.0,
    );
}

fn draw_reads(
    elements: &mut Vec<VisualElement>,
    sample: &SamplePlotData,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    let (placed_reads, lanes) = place_reads(&sample.reads, opts.squash);
    if lanes > 10 {
        for lane in (10..lanes).step_by(10) {
            let line_y = y + lane as f64 * opts.lane_height as f64;
            elements.push(VisualElement::Line {
                x1: geom.x0,
                y1: line_y,
                x2: geom.x0 + geom.width,
                y2: line_y,
                color: Rgb(228, 228, 228),
                width: 1,
            });
        }
    }
    for placed in placed_reads {
        draw_read(
            elements,
            placed.read,
            plot,
            opts,
            geom,
            y + placed.lane as f64 * (if opts.squash { 3.0 } else { 16.0 }),
        );
    }
}

fn draw_read(
    elements: &mut Vec<VisualElement>,
    read: &ReadModel,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    let read_color = match read.haplotype {
        Some(1) => opts.style.haplotype_1,
        Some(2) => opts.style.haplotype_2,
        _ => opts.style.read_forward,
    };
    let h = if opts.squash { 3.0 } else { 12.0 };
    let visible_match_start = read
        .segments
        .iter()
        .filter_map(|segment| match segment {
            ReadSegment::Match { ref_start, .. } | ReadSegment::Mismatch { ref_start, .. } => {
                Some(*ref_start)
            }
            _ => None,
        })
        .min();
    let visible_match_end = read
        .segments
        .iter()
        .filter_map(|segment| match segment {
            ReadSegment::Match { ref_start, len, .. }
            | ReadSegment::Mismatch { ref_start, len, .. } => Some(*ref_start + *len),
            _ => None,
        })
        .max();
    for segment in &read.segments {
        match segment {
            ReadSegment::Match { ref_start, len, .. } => push_read_body(
                elements,
                *ref_start,
                *ref_start + *len,
                y,
                h,
                read_color,
                read.is_reverse,
                visible_match_start,
                visible_match_end,
                plot,
                geom,
            ),
            ReadSegment::Mismatch {
                ref_start,
                len,
                base,
                ..
            } => push_read_rect(
                elements,
                *ref_start,
                *ref_start + *len,
                y,
                h,
                base_color(*base),
                plot,
                geom,
            ),
            ReadSegment::Del { ref_start, len } => elements.push(VisualElement::Line {
                x1: geom.x(*ref_start, plot),
                y1: y + h / 2.0,
                x2: geom.x(*ref_start + *len, plot),
                y2: y + h / 2.0,
                color: opts.style.deletion.into(),
                width: 2,
            }),
            ReadSegment::Ins { ref_pos, .. } => elements.push(VisualElement::Triangle {
                x: geom.x(*ref_pos, plot),
                y: y + h / 2.0,
                size: 5.0,
                color: opts.style.insertion.into(),
            }),
            ReadSegment::SoftClip { ref_pos, bases } => {
                let x1 = geom.x(*ref_pos, plot);
                let soft_w = (bases.len() as f64 * 1.5).clamp(2.0, 18.0);
                elements.push(VisualElement::Rect {
                    x1,
                    y1: y + h * 0.2,
                    x2: x1 + soft_w,
                    y2: y + h * 0.8,
                    color: lighten(read_color, 0.45).into(),
                });
            }
            ReadSegment::Skip { ref_start, len } => elements.push(VisualElement::Line {
                x1: geom.x(*ref_start, plot),
                y1: y + h / 2.0,
                x2: geom.x(*ref_start + *len, plot),
                y2: y + h / 2.0,
                color: if opts.squash {
                    Rgb(150, 150, 150)
                } else {
                    Rgb(205, 205, 205)
                },
                width: 1,
            }),
        }
    }
    for modification in &read.modifications {
        elements.push(VisualElement::Circle {
            x: geom.x(modification.ref_pos, plot),
            y: y - 1.0,
            radius: 2.0,
            color: opts.style.modification.into(),
        });
    }
}

fn push_read_rect(
    elements: &mut Vec<VisualElement>,
    start: i64,
    end: i64,
    y: f64,
    h: f64,
    color: RGBColor,
    plot: &RegionPlot,
    geom: PlotGeom,
) {
    let x1 = geom.x(start, plot);
    elements.push(VisualElement::Rect {
        x1,
        y1: y,
        x2: geom.x(end, plot).max(x1 + 1.0),
        y2: y + h,
        color: color.into(),
    });
}

fn push_read_body(
    elements: &mut Vec<VisualElement>,
    start: i64,
    end: i64,
    y: f64,
    h: f64,
    color: RGBColor,
    is_reverse: bool,
    visible_match_start: Option<i64>,
    visible_match_end: Option<i64>,
    plot: &RegionPlot,
    geom: PlotGeom,
) {
    let x1 = geom.x(start, plot);
    let x2 = geom.x(end, plot).max(x1 + 1.0);
    let arrow_px = 7.0_f64.min((x2 - x1) * 0.25);
    let y_mid = y + h / 2.0;
    let use_right_arrow = !is_reverse && visible_match_end == Some(end) && arrow_px >= 2.0;
    let use_left_arrow = is_reverse && visible_match_start == Some(start) && arrow_px >= 2.0;

    if use_right_arrow {
        if x2 - x1 > arrow_px {
            elements.push(VisualElement::Rect {
                x1,
                y1: y,
                x2: x2 - arrow_px,
                y2: y + h,
                color: color.into(),
            });
        }
        elements.push(VisualElement::Polygon {
            points: vec![(x2 - arrow_px, y), (x2, y_mid), (x2 - arrow_px, y + h)],
            color: color.into(),
        });
    } else if use_left_arrow {
        if x2 - x1 > arrow_px {
            elements.push(VisualElement::Rect {
                x1: x1 + arrow_px,
                y1: y,
                x2,
                y2: y + h,
                color: color.into(),
            });
        }
        elements.push(VisualElement::Polygon {
            points: vec![(x1 + arrow_px, y), (x1, y_mid), (x1 + arrow_px, y + h)],
            color: color.into(),
        });
    } else {
        elements.push(VisualElement::Rect {
            x1,
            y1: y,
            x2,
            y2: y + h,
            color: color.into(),
        });
    }
}

fn coverage_depths(sample: &SamplePlotData, region_start: i64) -> Vec<(i64, u32)> {
    if !sample.pileup.is_empty() {
        sample
            .pileup
            .iter()
            .enumerate()
            .map(|(idx, pile)| (region_start + idx as i64, pile.depth()))
            .collect()
    } else {
        sample
            .coverage
            .iter()
            .map(|point| (point.pos, point.depth))
            .collect()
    }
}

fn base_color(base: u8) -> RGBColor {
    match base.to_ascii_uppercase() {
        b'A' => RGBColor(0, 255, 0),
        b'C' => RGBColor(0, 0, 255),
        b'G' => RGBColor(235, 140, 0),
        b'T' => RGBColor(255, 0, 0),
        _ => RGBColor(130, 130, 130),
    }
}

fn lighten(color: RGBColor, amount: f64) -> RGBColor {
    let blend = |channel: u8| {
        (channel as f64 + (255.0 - channel as f64) * amount)
            .round()
            .clamp(0.0, 255.0) as u8
    };
    RGBColor(blend(color.0), blend(color.1), blend(color.2))
}
