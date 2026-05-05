use crate::layout::place_reads;
use crate::model::{GeneModel, ReadModel, ReadSegment, RegionPlot, SamplePlotData};
use crate::render::PlotOptions;
use plotters::style::RGBColor;

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
        y += opts.gene_height as f64;
    }

    if let Some(reference) = &plot.reference {
        draw_reference_strip(&mut elements, reference, plot, opts, geom, y);
        y += opts.reference_height as f64;
    }

    for sample in &plot.samples {
        draw_sample(&mut elements, sample, plot, opts, geom, y);
        y += sample_height(sample, opts) as f64;
    }

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
    let gene_h = if plot.genes.is_empty() {
        0
    } else {
        opts.gene_height
    };
    let reference_h = if plot.reference.is_some() {
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
    let (_, lanes) = place_reads(&sample.reads);
    opts.sample_label_height
        + opts.coverage_height
        + opts.base_track_height
        + opts.read_track_gap
        + lanes.max(1) as u32 * opts.lane_height
        + opts.sample_gap
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
        y: y + 22.0,
        text: "Genes".to_string(),
        color: opts.style.text.into(),
        size: 13,
    });
    for (idx, gene) in genes.iter().enumerate() {
        let lane_y = y + 18.0 + (idx % 2) as f64 * 18.0;
        let x1 = geom.x(gene.start, plot);
        let x2 = geom.x(gene.end, plot).max(x1 + 1.0);
        elements.push(VisualElement::Line {
            x1,
            y1: lane_y,
            x2,
            y2: lane_y,
            color: opts.style.gene.into(),
            width: 2,
        });
        for &(start, end) in &gene.exons {
            let ex1 = geom.x(start, plot);
            elements.push(VisualElement::Rect {
                x1: ex1,
                y1: lane_y - 5.0,
                x2: geom.x(end, plot).max(ex1 + 2.0),
                y2: lane_y + 5.0,
                color: opts.style.exon.into(),
            });
        }
        elements.push(VisualElement::Text {
            x: x1,
            y: lane_y - 8.0,
            text: gene.name.clone(),
            color: opts.style.text.into(),
            size: 11,
        });
    }
}

fn draw_reference_strip(
    elements: &mut Vec<VisualElement>,
    reference: &[u8],
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    elements.push(VisualElement::Text {
        x: 12.0,
        y: y + 18.0,
        text: "Reference".to_string(),
        color: opts.style.text.into(),
        size: 13,
    });
    let px_per_base = geom.width / plot.span() as f64;
    let draw_text = px_per_base > 8.0;
    for (idx, base) in reference.iter().enumerate() {
        let pos = plot.start + idx as i64;
        if pos >= plot.end {
            break;
        }
        if *base == b'N' {
            continue;
        }
        let x1 = geom.x(pos, plot);
        let x2 = geom.base_x2(pos, plot);
        elements.push(VisualElement::Rect {
            x1,
            y1: y + 4.0,
            x2,
            y2: y + opts.reference_height as f64 - 4.0,
            color: base_color(*base).into(),
        });
        if draw_text {
            elements.push(VisualElement::Text {
                x: x1 + (x2 - x1) * 0.32,
                y: y + opts.reference_height as f64 - 8.0,
                text: (*base as char).to_string(),
                color: Rgb(255, 255, 255),
                size: 12,
            });
        }
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
    draw_sample_base_track(elements, sample, plot, opts, geom, base_y);
    let read_y = base_y + opts.base_track_height as f64 + opts.read_track_gap as f64;
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
    sample: &SamplePlotData,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    if sample.pileup.is_empty() {
        return;
    }
    let reference = plot.reference.as_ref();
    let px_per_base = geom.width / plot.span() as f64;
    let draw_text = px_per_base > 10.0;
    for (idx, pile) in sample.pileup.iter().enumerate() {
        let pos = plot.start + idx as i64;
        if pos >= plot.end {
            break;
        }
        let ref_base = reference
            .and_then(|seq| seq.get(idx))
            .copied()
            .unwrap_or(b'N');
        let total = pile.depth().max(1) as f32;
        let mut lanes = Vec::new();
        for base in [b'A', b'C', b'G', b'T'] {
            let count = pile.base_count(base);
            if count > 0
                && (reference.is_none() || base == ref_base || (count as f32 / total) >= opts.min_alt_af)
            {
                lanes.push(base);
            }
        }
        if lanes.is_empty() && reference.is_some() {
            lanes.push(ref_base);
        } else if lanes.is_empty() {
            continue;
        }
        let lane_h = (opts.base_track_height as f64 - 4.0) / lanes.len().max(1) as f64;
        for (lane_idx, base) in lanes.into_iter().enumerate() {
            let x1 = geom.x(pos, plot);
            let x2 = geom.base_x2(pos, plot);
            let y1 = y + 2.0 + lane_idx as f64 * lane_h;
            elements.push(VisualElement::Rect {
                x1,
                y1,
                x2,
                y2: y1 + lane_h.max(1.0),
                color: base_color(base).into(),
            });
            if draw_text {
                elements.push(VisualElement::Text {
                    x: x1 + (x2 - x1) * 0.30,
                    y: y1 + lane_h - 2.0,
                    text: (base as char).to_string(),
                    color: Rgb(255, 255, 255),
                    size: 10,
                });
            }
        }
    }
}

fn draw_reads(
    elements: &mut Vec<VisualElement>,
    sample: &SamplePlotData,
    plot: &RegionPlot,
    opts: &PlotOptions,
    geom: PlotGeom,
    y: f64,
) {
    let (placed_reads, lanes) = place_reads(&sample.reads);
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
            y + placed.lane as f64 * opts.lane_height as f64,
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
    let h = (opts.lane_height as f64 * 0.62).max(3.0);
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
