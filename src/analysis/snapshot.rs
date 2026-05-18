use crate::analysis::bam_scan::find_bai_path;
use crate::io::annotation::AnnotationFormat;
use crate::io::text::open_maybe_gz;
use anyhow::{anyhow, bail, Context, Result};
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::{bam, fasta, sam};
use regex::Regex;
use region_plot::{
    render_to_path, BasePileup, CoveragePoint, GeneModel, MarkerType, PlotOptions, ReadModel,
    ReadSegment, RegionPlot, SamplePlotData, SnappingMarker,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::Path;

const DEFAULT_SINGLE_POS_FLANK_BP: u64 = 500;
const MAX_REGION_BP: u64 = 1_000_000;
const REFERENCE_DISPLAY_MAX_BP: u64 = 1_000_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicRegion {
    pub chrom: String,
    /// 0-based inclusive.
    pub start: u64,
    /// 0-based exclusive.
    pub end: u64,
}

impl GenomicRegion {
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    fn to_query_string(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start + 1, self.end)
    }
}

#[derive(Debug, Clone)]
pub struct SnapshotConfig {
    pub bam_path: String,
    pub bai_path: Option<String>,
    pub region: GenomicRegion,
    pub annotation_path: Option<String>,
    pub annotation_format: AnnotationFormat,
    pub reference_path: Option<String>,
    pub output_path: String,
    pub mapq_threshold: u8,
    pub max_reads: usize,
    pub width: u32,
    pub min_height: u32,
    pub show_reference: bool,
    pub show_genes: bool,
    pub show_reference_base_track: bool,
    pub show_sample_base_track: bool,
    pub squash: bool,
    pub markers_path: Option<String>,
    pub inline_markers: Vec<region_plot::SnappingMarker>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnapshotOutputFormat {
    Auto,
    Png,
    Svg,
}

pub fn run_snapshot(config: &SnapshotConfig) -> Result<()> {
    let plot = build_region_plot(config)?;
    let opts = PlotOptions {
        width: config.width,
        min_height: config.min_height,
        font_family: "sans-serif".to_string(),
        squash: config.squash,
        show_reference_base_track: config.show_reference_base_track,
        show_sample_base_track: config.show_sample_base_track,
        ..PlotOptions::default()
    };
    render_to_path(&plot, &opts, &config.output_path)
        .with_context(|| format!("failed to render {}", config.output_path))?;
    Ok(())
}

pub fn build_region_plot(config: &SnapshotConfig) -> Result<RegionPlot> {
    if config.region.len() == 0 {
        bail!("snapshot region has zero length");
    }
    if config.region.len() > MAX_REGION_BP {
        bail!(
            "snapshot region is {} bp; v1 supports regions up to {} bp",
            config.region.len(),
            MAX_REGION_BP
        );
    }

    let reference = if config.show_reference {
        load_reference(config)?
    } else {
        None
    };
    let (reads, coverage, pileup) = extract_bam_snapshot(config, reference.as_deref())?;
    let genes = if config.show_genes {
        load_genes(config)?
    } else {
        Vec::new()
    };
    let sample_name = Path::new(&config.bam_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(&config.bam_path)
        .to_string();

    let mut markers = config.inline_markers.clone();
    if let Some(ref path) = config.markers_path {
        let file_markers = parse_markers_file(path)?;
        markers.extend(file_markers);
    }

    Ok(RegionPlot {
        chrom: config.region.chrom.clone(),
        start: config.region.start as i64,
        end: config.region.end as i64,
        reference,
        genes,
        samples: vec![SamplePlotData {
            name: sample_name,
            reads,
            pileup,
            coverage,
        }],
        markers,
    })
}

pub fn parse_markers_file(path: &str) -> Result<Vec<SnappingMarker>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open markers file: {path}"))?;
    let reader = std::io::BufReader::new(file);
    let mut markers = Vec::new();

    for (line_idx, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() < 4 {
            bail!("invalid marker line {}: expected at least 4 tab-separated columns (chrom, position, label, marker_type)", line_idx + 1);
        }

        let pos: i64 = parts[1].parse()
            .with_context(|| format!("invalid marker position '{}' at line {}", parts[1], line_idx + 1))?;
        let label = parts[2].to_string();
        let marker_type_str = parts[3].to_lowercase();
        let marker_type = match marker_type_str.as_str() {
            "variant" | "snv" | "indel" => MarkerType::Variant,
            "structuralvariant" | "sv" => MarkerType::StructuralVariant,
            "regionofinterest" | "roi" | "region" => MarkerType::RegionOfInterest,
            other => bail!("unknown marker type '{}' at line {}. Supported: Variant, StructuralVariant, RegionOfInterest", other, line_idx + 1),
        };
        let end_pos = if parts.len() > 4 && !parts[4].is_empty() {
            let ep: i64 = parts[4].parse()
                .with_context(|| format!("invalid marker end_pos '{}' at line {}", parts[4], line_idx + 1))?;
            Some(ep)
        } else {
            None
        };

        markers.push(SnappingMarker {
            pos,
            label,
            marker_type,
            end_pos,
        });
    }

    Ok(markers)
}

pub fn resolve_snapshot_region(
    raw_region: &str,
    annotation_path: Option<&str>,
    annotation_format: AnnotationFormat,
) -> Result<GenomicRegion> {
    resolve_snapshot_region_with_flank(raw_region, annotation_path, annotation_format, None)
}

pub fn resolve_snapshot_region_with_flank(
    raw_region: &str,
    annotation_path: Option<&str>,
    annotation_format: AnnotationFormat,
    flank_bp: Option<u64>,
) -> Result<GenomicRegion> {
    if let Ok(region) = parse_region(raw_region) {
        return Ok(region);
    }

    let Some(path) = annotation_path else {
        bail!("region must look like chr:start-end or a gene name when --annotation is provided");
    };

    lookup_gene_region(path, raw_region, annotation_format)
        .map(|coords| pad_gene_region(coords, flank_bp))
        .ok_or_else(|| anyhow!("invalid region or gene not found: {raw_region}"))
}

pub fn parse_region(raw: &str) -> Result<GenomicRegion> {
    let (chrom, rest) = raw
        .split_once(':')
        .ok_or_else(|| anyhow!("region must look like chr:start-end or chr:pos"))?;
    if chrom.is_empty() {
        bail!("region chromosome is empty");
    }

    let rest = rest.replace(',', "");
    if let Some((start_raw, end_raw)) = rest.split_once('-') {
        let start_1: u64 = start_raw
            .parse()
            .with_context(|| format!("invalid region start: {start_raw}"))?;
        let end_1: u64 = end_raw
            .parse()
            .with_context(|| format!("invalid region end: {end_raw}"))?;
        if start_1 == 0 || end_1 < start_1 {
            bail!("invalid region coordinates: {raw}");
        }
        Ok(GenomicRegion {
            chrom: chrom.to_string(),
            start: start_1 - 1,
            end: end_1,
        })
    } else {
        let pos_1: u64 = rest
            .parse()
            .with_context(|| format!("invalid region position: {rest}"))?;
        if pos_1 == 0 {
            bail!("region positions are 1-based and must be positive");
        }
        let center_0 = pos_1 - 1;
        Ok(GenomicRegion {
            chrom: chrom.to_string(),
            start: center_0.saturating_sub(DEFAULT_SINGLE_POS_FLANK_BP),
            end: center_0 + DEFAULT_SINGLE_POS_FLANK_BP + 1,
        })
    }
}

fn pad_gene_region(
    (chrom, start, end): (String, i64, i64),
    flank_bp: Option<u64>,
) -> GenomicRegion {
    let padding = flank_bp.map(|v| v as i64).unwrap_or_else(|| {
        let span = (end - start).abs().max(1);
        (span / 8).clamp(200, 5_000)
    });
    GenomicRegion {
        chrom,
        start: (start - padding).max(1) as u64,
        end: (end + padding).max(start + 1) as u64,
    }
}

fn lookup_gene_region(
    path: &str,
    gene_name: &str,
    annotation_format: AnnotationFormat,
) -> Option<(String, i64, i64)> {
    let gene_name_lower = gene_name.trim().to_lowercase();
    if gene_name_lower.is_empty() {
        return None;
    }

    match annotation_format {
        AnnotationFormat::Auto => {
            if path.ends_with(".gtf") || path.ends_with(".gtf.gz") {
                lookup_gtf_gene_region(path, &gene_name_lower)
                    .or_else(|| lookup_bed_gene_region(path, &gene_name_lower))
            } else {
                lookup_bed_gene_region(path, &gene_name_lower)
                    .or_else(|| lookup_gtf_gene_region(path, &gene_name_lower))
            }
        }
        AnnotationFormat::Gtf => lookup_gtf_gene_region(path, &gene_name_lower),
        AnnotationFormat::Bed12 => lookup_bed_gene_region(path, &gene_name_lower),
    }
}

fn lookup_gtf_gene_region(path: &str, gene_name_lower: &str) -> Option<(String, i64, i64)> {
    let reader = open_maybe_gz(path).ok()?;
    let gene_re = Regex::new(r#"gene_name "([^"]+)""#).ok()?;
    for line in reader.lines() {
        let line = line.ok()?;
        if line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 || parts[2] != "gene" {
            continue;
        }
        let caps = gene_re.captures(parts[8])?;
        let name = caps.get(1)?.as_str();
        if name.to_lowercase() != gene_name_lower {
            continue;
        }
        let start = parts[3].parse::<i64>().ok()?;
        let end = parts[4].parse::<i64>().ok()?;
        return Some((parts[0].to_string(), start, end));
    }
    None
}

fn lookup_bed_gene_region(path: &str, gene_name_lower: &str) -> Option<(String, i64, i64)> {
    let reader = open_maybe_gz(path).ok()?;
    for line in reader.lines() {
        let line = line.ok()?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let name = parts.get(3).copied().unwrap_or("").trim().to_lowercase();
        if name != gene_name_lower {
            continue;
        }
        let start = parts[1].parse::<i64>().ok()?;
        let end = parts[2].parse::<i64>().ok()?;
        return Some((parts[0].to_string(), start, end));
    }
    None
}

pub fn validate_output_format(path: &str, format: SnapshotOutputFormat) -> Result<()> {
    match format {
        SnapshotOutputFormat::Auto => match Path::new(path).extension().and_then(|e| e.to_str()) {
            Some(ext) if ext.eq_ignore_ascii_case("png") || ext.eq_ignore_ascii_case("svg") => {
                Ok(())
            }
            _ => bail!("--format auto requires an output path ending in .png or .svg"),
        },
        SnapshotOutputFormat::Png => require_ext(path, "png"),
        SnapshotOutputFormat::Svg => require_ext(path, "svg"),
    }
}

fn require_ext(path: &str, expected: &str) -> Result<()> {
    match Path::new(path).extension().and_then(|e| e.to_str()) {
        Some(ext) if ext.eq_ignore_ascii_case(expected) => Ok(()),
        _ => bail!("--format {expected} requires an output path ending in .{expected}"),
    }
}

fn extract_bam_snapshot(
    config: &SnapshotConfig,
    reference: Option<&[u8]>,
) -> Result<(Vec<ReadModel>, Vec<CoveragePoint>, Vec<BasePileup>)> {
    let bai_path = config
        .bai_path
        .clone()
        .or_else(|| find_bai_path(&config.bam_path))
        .ok_or_else(|| anyhow!("could not find BAI for {}", config.bam_path))?;
    let bai = bam::bai::read(&bai_path)
        .with_context(|| format!("failed to read BAI index {bai_path}"))?;
    let file = File::open(&config.bam_path)
        .with_context(|| format!("could not open BAM {}", config.bam_path))?;
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(bai)
        .build_from_reader(file)?;
    let header = reader.read_header()?;
    let fetch_chrom = chrom_aliases(&config.region.chrom)
        .into_iter()
        .find(|alias| header.reference_sequences().get(alias.as_bytes()).is_some())
        .ok_or_else(|| anyhow!("chromosome {} not found in BAM header", config.region.chrom))?;

    let region_str = format!(
        "{}:{}-{}",
        fetch_chrom,
        config.region.start + 1,
        config.region.end
    );
    let query_region: noodles::core::Region = region_str.parse()?;
    let query = reader.query(&header, &query_region)?;

    let mut displayed_reads = Vec::new();
    let mut coverage = vec![0_u32; config.region.len() as usize];
    let mut pileup = vec![BasePileup::default(); config.region.len() as usize];
    let mut seen = 0_usize;
    let mut rng = 0x9e37_79b9_7f4a_7c15_u64;

    for result in query {
        let record = result?;
        if !record_passes(&record, config.mapq_threshold) {
            continue;
        }
        add_record_coverage(&record, &config.region, &mut coverage);
        add_record_pileup(&record, &config.region, &mut pileup);
        if let Some(read) = record_to_read_model(&header, &record, &config.region, reference) {
            seen += 1;
            if displayed_reads.len() < config.max_reads {
                displayed_reads.push(read);
            } else if config.max_reads > 0 {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                let idx = (rng as usize) % seen;
                if idx < config.max_reads {
                    displayed_reads[idx] = read;
                }
            }
        }
    }

    if seen > displayed_reads.len() {
        eprintln!(
            "warning: displayed {} of {} passing reads in {}",
            displayed_reads.len(),
            seen,
            config.region.to_query_string()
        );
    }

    let coverage = coverage
        .into_iter()
        .enumerate()
        .map(|(i, depth)| CoveragePoint {
            pos: config.region.start as i64 + i as i64,
            depth,
        })
        .collect();

    Ok((displayed_reads, coverage, pileup))
}

fn record_passes(record: &bam::Record, mapq_threshold: u8) -> bool {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        return false;
    }
    let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
    mapq >= mapq_threshold
}

pub fn record_to_read_model(
    _header: &sam::Header,
    record: &bam::Record,
    region: &GenomicRegion,
    reference: Option<&[u8]>,
) -> Option<ReadModel> {
    let start = record.alignment_start()?.ok()?.get() as u64 - 1;
    let bases: Vec<u8> = record.sequence().iter().collect();
    let qualities = record.quality_scores().as_ref().to_vec();
    let segments = cigar_to_segments(
        start,
        record.cigar().iter().filter_map(Result::ok).map(|op| {
            let kind = match op.kind() {
                Kind::Match => SnapshotCigarKind::Match,
                Kind::SequenceMatch => SnapshotCigarKind::Match,
                Kind::SequenceMismatch => SnapshotCigarKind::Match,
                Kind::Insertion => SnapshotCigarKind::Ins,
                Kind::Deletion => SnapshotCigarKind::Del,
                Kind::Skip => SnapshotCigarKind::Skip,
                Kind::SoftClip => SnapshotCigarKind::SoftClip,
                Kind::HardClip => SnapshotCigarKind::HardClip,
                Kind::Pad => SnapshotCigarKind::Pad,
            };
            SnapshotCigarOp {
                kind,
                len: op.len(),
            }
        }),
        &bases,
        region,
        reference,
    );

    let mut ref_end = start;
    for op in record.cigar().iter().filter_map(Result::ok) {
        if matches!(
            op.kind(),
            Kind::Match
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Skip
        ) {
            ref_end += op.len() as u64;
        }
    }

    Some(ReadModel {
        name: record
            .name()
            .map(|name| String::from_utf8_lossy(name.as_ref()).to_string())
            .unwrap_or_else(|| "*".to_string()),
        start: start as i64,
        end: ref_end.max(start + 1) as i64,
        is_reverse: record.flags().is_reverse_complemented(),
        mapq: record.mapping_quality().map(|m| m.get()).unwrap_or(255),
        segments,
        bases: Some(bases),
        qualities: Some(qualities),
        haplotype: hp_tag(record),
        modifications: Vec::new(), // TODO: convert MM/ML base modification tags.
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SnapshotCigarKind {
    Match,
    Ins,
    Del,
    Skip,
    SoftClip,
    HardClip,
    Pad,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SnapshotCigarOp {
    pub kind: SnapshotCigarKind,
    pub len: usize,
}

pub fn cigar_to_segments<I>(
    ref_start: u64,
    ops: I,
    bases: &[u8],
    region: &GenomicRegion,
    reference: Option<&[u8]>,
) -> Vec<ReadSegment>
where
    I: IntoIterator<Item = SnapshotCigarOp>,
{
    let mut ref_pos = ref_start;
    let mut query_pos = 0_usize;
    let mut segments = Vec::new();

    for op in ops {
        match op.kind {
            SnapshotCigarKind::Match => {
                let op_start = ref_pos;
                let op_end = ref_pos + op.len as u64;
                let clip_start = op_start.max(region.start);
                let clip_end = op_end.min(region.end);
                if clip_end > clip_start {
                    push_match_or_mismatch_segments(
                        &mut segments,
                        clip_start,
                        clip_end,
                        op_start,
                        query_pos,
                        bases,
                        region,
                        reference,
                    );
                }
                ref_pos = op_end;
                query_pos += op.len;
            }
            SnapshotCigarKind::Ins => {
                if region.start <= ref_pos && ref_pos <= region.end {
                    segments.push(ReadSegment::Ins {
                        ref_pos: ref_pos as i64,
                        bases: slice_bases(bases, query_pos, op.len),
                    });
                }
                query_pos += op.len;
            }
            SnapshotCigarKind::Del => {
                let op_start = ref_pos;
                let op_end = ref_pos + op.len as u64;
                let clip_start = op_start.max(region.start);
                let clip_end = op_end.min(region.end);
                if clip_end > clip_start {
                    segments.push(ReadSegment::Del {
                        ref_start: clip_start as i64,
                        len: (clip_end - clip_start) as i64,
                    });
                }
                ref_pos = op_end;
            }
            SnapshotCigarKind::Skip => {
                let op_start = ref_pos;
                let op_end = ref_pos + op.len as u64;
                let clip_start = op_start.max(region.start);
                let clip_end = op_end.min(region.end);
                if clip_end > clip_start {
                    segments.push(ReadSegment::Skip {
                        ref_start: clip_start as i64,
                        len: (clip_end - clip_start) as i64,
                    });
                }
                ref_pos = op_end;
            }
            SnapshotCigarKind::SoftClip => {
                if region.start <= ref_pos && ref_pos <= region.end {
                    segments.push(ReadSegment::SoftClip {
                        ref_pos: ref_pos as i64,
                        bases: slice_bases(bases, query_pos, op.len),
                    });
                }
                query_pos += op.len;
            }
            SnapshotCigarKind::HardClip | SnapshotCigarKind::Pad => {}
        }
    }

    segments
}

fn push_match_or_mismatch_segments(
    segments: &mut Vec<ReadSegment>,
    clip_start: u64,
    clip_end: u64,
    op_start: u64,
    query_pos: usize,
    bases: &[u8],
    region: &GenomicRegion,
    reference: Option<&[u8]>,
) {
    let Some(reference) = reference else {
        segments.push(ReadSegment::Match {
            ref_start: clip_start as i64,
            len: (clip_end - clip_start) as i64,
            query_start: query_pos + (clip_start - op_start) as usize,
        });
        return;
    };

    let mut run_start = clip_start;
    let mut run_query_start = query_pos + (clip_start - op_start) as usize;
    let mut run_is_match = true;
    let mut run_base = b'N';

    for pos in clip_start..clip_end {
        let q_idx = query_pos + (pos - op_start) as usize;
        let read_base = bases
            .get(q_idx)
            .copied()
            .unwrap_or(b'N')
            .to_ascii_uppercase();
        let ref_idx = (pos - region.start) as usize;
        let ref_base = reference
            .get(ref_idx)
            .copied()
            .unwrap_or(b'N')
            .to_ascii_uppercase();
        let is_match = read_base == ref_base || read_base == b'N' || ref_base == b'N';

        if pos == clip_start {
            run_is_match = is_match;
            run_base = read_base;
        } else if is_match != run_is_match || (!is_match && read_base != run_base) {
            push_alignment_run(
                segments,
                run_start,
                pos,
                run_query_start,
                run_is_match,
                run_base,
            );
            run_start = pos;
            run_query_start = q_idx;
            run_is_match = is_match;
            run_base = read_base;
        }
    }

    push_alignment_run(
        segments,
        run_start,
        clip_end,
        run_query_start,
        run_is_match,
        run_base,
    );
}

fn push_alignment_run(
    segments: &mut Vec<ReadSegment>,
    start: u64,
    end: u64,
    query_start: usize,
    is_match: bool,
    base: u8,
) {
    if end <= start {
        return;
    }
    if is_match {
        segments.push(ReadSegment::Match {
            ref_start: start as i64,
            len: (end - start) as i64,
            query_start,
        });
    } else {
        segments.push(ReadSegment::Mismatch {
            ref_start: start as i64,
            len: (end - start) as i64,
            query_start,
            base,
        });
    }
}

fn slice_bases(bases: &[u8], start: usize, len: usize) -> Vec<u8> {
    bases
        .get(start..start.saturating_add(len).min(bases.len()))
        .unwrap_or_default()
        .to_vec()
}

fn hp_tag(record: &bam::Record) -> Option<u8> {
    let tag = Tag::new(b'H', b'P');
    match record.data().get(&tag)?.ok()? {
        Value::UInt8(v) => Some(v),
        Value::Int8(v) => u8::try_from(v).ok(),
        Value::UInt16(v) => u8::try_from(v).ok(),
        Value::Int16(v) => u8::try_from(v).ok(),
        Value::UInt32(v) => u8::try_from(v).ok(),
        Value::Int32(v) => u8::try_from(v).ok(),
        _ => None,
    }
}

pub fn add_coverage_from_ops<I>(
    ref_start: u64,
    ops: I,
    region: &GenomicRegion,
    coverage: &mut [u32],
) where
    I: IntoIterator<Item = SnapshotCigarOp>,
{
    let mut ref_pos = ref_start;
    for op in ops {
        match op.kind {
            SnapshotCigarKind::Match => {
                let op_start = ref_pos;
                let op_end = ref_pos + op.len as u64;
                let start = op_start.max(region.start);
                let end = op_end.min(region.end);
                if end > start {
                    for idx in (start - region.start) as usize..(end - region.start) as usize {
                        coverage[idx] += 1;
                    }
                }
                ref_pos = op_end;
            }
            SnapshotCigarKind::Del | SnapshotCigarKind::Skip => {
                ref_pos += op.len as u64;
            }
            SnapshotCigarKind::Ins | SnapshotCigarKind::SoftClip => {}
            SnapshotCigarKind::HardClip | SnapshotCigarKind::Pad => {}
        }
    }
}

fn add_record_coverage(record: &bam::Record, region: &GenomicRegion, coverage: &mut [u32]) {
    let Some(start) = record
        .alignment_start()
        .and_then(|p| p.ok())
        .map(|p| p.get() as u64 - 1)
    else {
        return;
    };
    add_coverage_from_ops(
        start,
        record
            .cigar()
            .iter()
            .filter_map(Result::ok)
            .map(|op| SnapshotCigarOp {
                kind: match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        SnapshotCigarKind::Match
                    }
                    Kind::Insertion => SnapshotCigarKind::Ins,
                    Kind::Deletion => SnapshotCigarKind::Del,
                    Kind::Skip => SnapshotCigarKind::Skip,
                    Kind::SoftClip => SnapshotCigarKind::SoftClip,
                    Kind::HardClip => SnapshotCigarKind::HardClip,
                    Kind::Pad => SnapshotCigarKind::Pad,
                },
                len: op.len(),
            }),
        region,
        coverage,
    );
}

fn add_record_pileup(record: &bam::Record, region: &GenomicRegion, pileup: &mut [BasePileup]) {
    let Some(start) = record
        .alignment_start()
        .and_then(|p| p.ok())
        .map(|p| p.get() as u64 - 1)
    else {
        return;
    };
    let bases: Vec<u8> = record.sequence().iter().collect();
    add_pileup_from_ops(
        start,
        record
            .cigar()
            .iter()
            .filter_map(Result::ok)
            .map(|op| SnapshotCigarOp {
                kind: match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        SnapshotCigarKind::Match
                    }
                    Kind::Insertion => SnapshotCigarKind::Ins,
                    Kind::Deletion => SnapshotCigarKind::Del,
                    Kind::Skip => SnapshotCigarKind::Skip,
                    Kind::SoftClip => SnapshotCigarKind::SoftClip,
                    Kind::HardClip => SnapshotCigarKind::HardClip,
                    Kind::Pad => SnapshotCigarKind::Pad,
                },
                len: op.len(),
            }),
        &bases,
        region,
        pileup,
    );
}

pub fn add_pileup_from_ops<I>(
    ref_start: u64,
    ops: I,
    bases: &[u8],
    region: &GenomicRegion,
    pileup: &mut [BasePileup],
) where
    I: IntoIterator<Item = SnapshotCigarOp>,
{
    let mut ref_pos = ref_start;
    let mut query_pos = 0_usize;

    for op in ops {
        match op.kind {
            SnapshotCigarKind::Match => {
                let op_start = ref_pos;
                for offset in 0..op.len {
                    let pos = op_start + offset as u64;
                    if pos >= region.start && pos < region.end {
                        let idx = (pos - region.start) as usize;
                        let base = bases.get(query_pos + offset).copied().unwrap_or(b'N');
                        if let Some(slot) = pileup.get_mut(idx) {
                            add_base_to_pileup(slot, base);
                        }
                    }
                }
                ref_pos += op.len as u64;
                query_pos += op.len;
            }
            SnapshotCigarKind::Ins => {
                if ref_pos >= region.start && ref_pos < region.end {
                    let idx = (ref_pos - region.start) as usize;
                    if let Some(slot) = pileup.get_mut(idx) {
                        slot.ins += 1;
                        slot.ins_len += op.len as u32;
                    }
                }
                query_pos += op.len;
            }
            SnapshotCigarKind::Del => {
                let op_start = ref_pos;
                for offset in 0..op.len {
                    let pos = op_start + offset as u64;
                    if pos >= region.start && pos < region.end {
                        let idx = (pos - region.start) as usize;
                        if let Some(slot) = pileup.get_mut(idx) {
                            slot.del += 1;
                            slot.del_spanning += 1;
                            if offset == 0 {
                                slot.del_starts += 1;
                            }
                        }
                    }
                }
                ref_pos += op.len as u64;
            }
            SnapshotCigarKind::Skip => {
                ref_pos += op.len as u64;
            }
            SnapshotCigarKind::SoftClip => {
                query_pos += op.len;
            }
            SnapshotCigarKind::HardClip | SnapshotCigarKind::Pad => {}
        }
    }
}

fn add_base_to_pileup(pile: &mut BasePileup, base: u8) {
    pile.total += 1;
    match base.to_ascii_uppercase() {
        b'A' => pile.a += 1,
        b'C' => pile.c += 1,
        b'G' => pile.g += 1,
        b'T' => pile.t += 1,
        _ => pile.n += 1,
    }
}

fn load_reference(config: &SnapshotConfig) -> Result<Option<Vec<u8>>> {
    let Some(path) = &config.reference_path else {
        return Ok(None);
    };
    if config.region.len() > REFERENCE_DISPLAY_MAX_BP {
        return Ok(None);
    }

    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("could not open indexed FASTA {path}"))?;

    let mut sequence = None;
    for alias in chrom_aliases(&config.region.chrom) {
        let region_str = format!(
            "{}:{}-{}",
            alias,
            config.region.start + 1,
            config.region.end
        );
        if let Ok(region) = region_str.parse() {
            if let Ok(rec) = reader.query(&region) {
                sequence = Some(rec.sequence().as_ref().to_vec());
                break;
            }
        }
    }

    Ok(sequence)
}

fn load_genes(config: &SnapshotConfig) -> Result<Vec<GeneModel>> {
    let Some(path) = &config.annotation_path else {
        return Ok(Vec::new());
    };
    let mut genes = Vec::new();
    for alias in chrom_aliases(&config.region.chrom) {
        let mut local_region = config.region.clone();
        local_region.chrom = alias;
        match config.annotation_format {
            AnnotationFormat::Auto => {
                if path.ends_with(".gtf") || path.ends_with(".gtf.gz") {
                    if let Ok(g) = load_gtf_genes(path, &local_region) {
                        if !g.is_empty() {
                            genes = g;
                            break;
                        }
                    }
                } else {
                    if let Ok(g) = load_bed12_genes(path, &local_region) {
                        if !g.is_empty() {
                            genes = g;
                            break;
                        }
                    }
                }
            }
            AnnotationFormat::Gtf => {
                if let Ok(g) = load_gtf_genes(path, &local_region) {
                    if !g.is_empty() {
                        genes = g;
                        break;
                    }
                }
            }
            AnnotationFormat::Bed12 => {
                if let Ok(g) = load_bed12_genes(path, &local_region) {
                    if !g.is_empty() {
                        genes = g;
                        break;
                    }
                }
            }
        }
    }
    Ok(genes)
}

#[derive(Debug, Default)]
struct GeneBuilder {
    name: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: Option<char>,
    exons: Vec<(u64, u64)>,
}

fn load_gtf_genes(path: &str, region: &GenomicRegion) -> Result<Vec<GeneModel>> {
    let reader = open_maybe_gz(path)?;
    let mut genes: HashMap<String, GeneBuilder> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 || fields[2] != "exon" || !chrom_match(fields[0], &region.chrom) {
            continue;
        }
        let start = fields[3].parse::<u64>()?.saturating_sub(1);
        let end = fields[4].parse::<u64>()?;
        if end <= region.start || start >= region.end {
            continue;
        }
        let attrs = parse_gtf_attrs(fields[8]);
        let gene_name = attrs
            .get("gene_name")
            .or_else(|| attrs.get("gene_id"))
            .cloned()
            .unwrap_or_else(|| "gene".to_string());
        let transcript_id = attrs
            .get("transcript_id")
            .or_else(|| attrs.get("gene_id"))
            .cloned()
            .unwrap_or_else(|| "unknown".to_string());

        let entry = genes
            .entry(transcript_id.clone())
            .or_insert_with(|| GeneBuilder {
                name: gene_name.clone(),
                chrom: fields[0].to_string(),
                start,
                end,
                strand: fields[6].chars().next(),
                exons: Vec::new(),
            });
        entry.start = entry.start.min(start);
        entry.end = entry.end.max(end);
        entry.exons.push((start, end));
    }

    Ok(build_gene_models(genes, region))
}

fn parse_gtf_attrs(raw: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    for part in raw.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let mut pieces = part.splitn(2, char::is_whitespace);
        if let (Some(key), Some(value)) = (pieces.next(), pieces.next()) {
            attrs.insert(key.to_string(), value.trim().trim_matches('"').to_string());
        }
    }
    attrs
}

fn load_bed12_genes(path: &str, region: &GenomicRegion) -> Result<Vec<GeneModel>> {
    let reader = open_maybe_gz(path)?;
    let mut genes = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 || !chrom_match(fields[0], &region.chrom) {
            continue;
        }
        let chrom_start: u64 = fields[1].parse()?;
        let chrom_end: u64 = fields[2].parse()?;
        if chrom_end <= region.start || chrom_start >= region.end {
            continue;
        }
        let name = fields.get(3).copied().unwrap_or("feature").to_string();
        let strand = fields.get(5).and_then(|s| s.chars().next());
        let mut exons = Vec::new();
        if fields.len() >= 12 {
            let sizes: Vec<u64> = fields[10]
                .trim_end_matches(',')
                .split(',')
                .filter_map(|s| s.parse().ok())
                .collect();
            let starts: Vec<u64> = fields[11]
                .trim_end_matches(',')
                .split(',')
                .filter_map(|s| s.parse().ok())
                .collect();
            for (size, rel_start) in sizes.into_iter().zip(starts) {
                let start = chrom_start + rel_start;
                let end = start + size;
                if end > region.start && start < region.end {
                    exons.push((start.max(region.start), end.min(region.end)));
                }
            }
        } else {
            exons.push((chrom_start.max(region.start), chrom_end.min(region.end)));
        }
        genes.insert(
            name.clone(),
            GeneBuilder {
                name,
                chrom: fields[0].to_string(),
                start: chrom_start,
                end: chrom_end,
                strand,
                exons,
            },
        );
    }

    Ok(build_gene_models(genes, region))
}

fn build_gene_models(
    genes: HashMap<String, GeneBuilder>,
    region: &GenomicRegion,
) -> Vec<GeneModel> {
    let mut out: Vec<_> = genes
        .into_values()
        .filter(|g| {
            chrom_match(&g.chrom, &region.chrom) && g.end > region.start && g.start < region.end
        })
        .map(|mut g| {
            g.exons.sort_unstable();
            g.exons.dedup();
            GeneModel {
                name: g.name,
                start: g.start.max(region.start) as i64,
                end: g.end.min(region.end) as i64,
                strand: g.strand,
                exons: g
                    .exons
                    .into_iter()
                    .map(|(s, e)| (s as i64, e as i64))
                    .collect(),
            }
        })
        .collect();
    out.sort_by_key(|g| (g.start, g.end, g.name.clone()));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use region_plot::{render_to_path, SamplePlotData};
    use tempfile::tempdir;

    fn region() -> GenomicRegion {
        GenomicRegion {
            chrom: "chr1".to_string(),
            start: 99,
            end: 130,
        }
    }

    #[test]
    fn parses_regions() {
        assert_eq!(
            parse_region("chr1:100-200").unwrap(),
            GenomicRegion {
                chrom: "chr1".to_string(),
                start: 99,
                end: 200
            }
        );
        assert_eq!(parse_region("chr1:100,000-101,000").unwrap().start, 99_999);
        let single = parse_region("chr1:100000").unwrap();
        assert_eq!(single.start, 99_499);
        assert_eq!(single.end, 100_500);
    }

    #[test]
    fn cigar_segments_cover_events() {
        let bases = b"AAAAACCGGGGGTTTTT";
        let r = region();
        assert_eq!(
            cigar_to_segments(
                100,
                [SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 10
                }],
                bases,
                &r,
                None
            )
            .len(),
            1
        );
        assert!(matches!(
            cigar_to_segments(
                100,
                [
                    SnapshotCigarOp {
                        kind: SnapshotCigarKind::Match,
                        len: 5
                    },
                    SnapshotCigarOp {
                        kind: SnapshotCigarKind::Ins,
                        len: 2
                    },
                    SnapshotCigarOp {
                        kind: SnapshotCigarKind::Match,
                        len: 5
                    },
                ],
                bases,
                &r,
                None
            )[1],
            ReadSegment::Ins { .. }
        ));
        assert!(matches!(
            cigar_to_segments(
                100,
                [
                    SnapshotCigarOp {
                        kind: SnapshotCigarKind::Match,
                        len: 5
                    },
                    SnapshotCigarOp {
                        kind: SnapshotCigarKind::Del,
                        len: 3
                    },
                    SnapshotCigarOp {
                        kind: SnapshotCigarKind::Match,
                        len: 5
                    },
                ],
                bases,
                &r,
                None
            )[1],
            ReadSegment::Del { .. }
        ));
        let soft = cigar_to_segments(
            100,
            [
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::SoftClip,
                    len: 5,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 10,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::SoftClip,
                    len: 5,
                },
            ],
            bases,
            &r,
            None,
        );
        assert!(soft
            .iter()
            .any(|s| matches!(s, ReadSegment::SoftClip { .. })));
        let skipped = cigar_to_segments(
            100,
            [
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 10,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Skip,
                    len: 100,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 10,
                },
            ],
            bases,
            &GenomicRegion {
                chrom: "chr1".to_string(),
                start: 100,
                end: 230,
            },
            None,
        );
        assert_eq!(skipped.len(), 3);
        assert!(skipped
            .iter()
            .any(|segment| matches!(segment, ReadSegment::Skip { .. })));
    }

    #[test]
    fn coverage_counts_only_matching_reference_bases() {
        let r = GenomicRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 115,
        };
        let mut cov = vec![0; r.len() as usize];
        add_coverage_from_ops(
            100,
            [
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 5,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Ins,
                    len: 2,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Del,
                    len: 3,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Skip,
                    len: 2,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 5,
                },
            ],
            &r,
            &mut cov,
        );
        assert_eq!(&cov[0..5], &[1, 1, 1, 1, 1]);
        assert_eq!(&cov[5..10], &[0, 0, 0, 0, 0]);
        assert_eq!(&cov[10..15], &[1, 1, 1, 1, 1]);
    }

    #[test]
    fn cigar_segments_split_mismatches_against_reference() {
        let r = GenomicRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 106,
        };
        let segments = cigar_to_segments(
            100,
            [SnapshotCigarOp {
                kind: SnapshotCigarKind::Match,
                len: 6,
            }],
            b"ACGTTC",
            &r,
            Some(b"ACGTAC"),
        );

        assert!(segments.iter().any(|s| matches!(
            s,
            ReadSegment::Mismatch {
                ref_start: 104,
                base: b'T',
                ..
            }
        )));
    }

    #[test]
    fn pileup_tracks_bases_and_cigar_events() {
        let r = GenomicRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 112,
        };
        let mut pileup = vec![BasePileup::default(); r.len() as usize];
        add_pileup_from_ops(
            100,
            [
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 4,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Ins,
                    len: 2,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 2,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Del,
                    len: 3,
                },
                SnapshotCigarOp {
                    kind: SnapshotCigarKind::Match,
                    len: 3,
                },
            ],
            b"ACGTGGACAAA",
            &r,
            &mut pileup,
        );

        assert_eq!(pileup[0].a, 1);
        assert_eq!(pileup[1].c, 1);
        assert_eq!(pileup[2].g, 1);
        assert_eq!(pileup[3].t, 1);
        assert_eq!(pileup[4].ins, 1);
        assert_eq!(pileup[4].ins_len, 2);
        assert_eq!(pileup[4].a, 1);
        assert_eq!(pileup[5].c, 1);
        assert_eq!(pileup[6].del_starts, 1);
        assert_eq!(pileup[6].del_spanning, 1);
        assert_eq!(pileup[7].del_spanning, 1);
        assert_eq!(pileup[8].del_spanning, 1);
        assert_eq!(pileup[9].a, 1);
    }

    #[test]
    fn renders_svg_smoke_test() {
        let tmp = tempdir().unwrap();
        let path = tmp.path().join("snap.svg");
        let plot = RegionPlot {
            chrom: "chr1".to_string(),
            start: 100,
            end: 120,
            reference: Some(b"ACGTACGTACGTACGTACGT".to_vec()),
            genes: vec![GeneModel {
                name: "GENE1".to_string(),
                start: 102,
                end: 118,
                strand: Some('+'),
                exons: vec![(102, 108), (112, 118)],
            }],
            samples: vec![SamplePlotData {
                name: "sample".to_string(),
                reads: Vec::new(),
                pileup: Vec::new(),
                coverage: (100..120)
                    .map(|pos| CoveragePoint { pos, depth: 1 })
                    .collect(),
            }],
            markers: Vec::new(),
        };
        render_to_path(&plot, &PlotOptions::default(), &path).unwrap();
        let svg = std::fs::read_to_string(path).unwrap();
        assert!(svg.contains("<svg"));
    }

    #[test]
    fn test_parse_markers_file() {
        let tmp = tempdir().unwrap();
        let path = tmp.path().join("markers.tsv");
        std::fs::write(&path, 
            "chr1\t10020\tC>T\tVariant\n\
             chr1\t10100\tDEL_1\tStructuralVariant\t10250\n\
             chr1\t10050\tROI_A\tRegionOfInterest\t10080\n"
        ).unwrap();

        let parsed = parse_markers_file(path.to_str().unwrap()).unwrap();
        assert_eq!(parsed.len(), 3);

        assert_eq!(parsed[0].pos, 10020);
        assert_eq!(parsed[0].label, "C>T");
        assert_eq!(parsed[0].marker_type, MarkerType::Variant);
        assert_eq!(parsed[0].end_pos, None);

        assert_eq!(parsed[1].pos, 10100);
        assert_eq!(parsed[1].label, "DEL_1");
        assert_eq!(parsed[1].marker_type, MarkerType::StructuralVariant);
        assert_eq!(parsed[1].end_pos, Some(10250));

        assert_eq!(parsed[2].pos, 10050);
        assert_eq!(parsed[2].label, "ROI_A");
        assert_eq!(parsed[2].marker_type, MarkerType::RegionOfInterest);
        assert_eq!(parsed[2].end_pos, Some(10080));
    }

    #[test]
    fn test_renders_svg_with_markers() {
        let tmp = tempdir().unwrap();
        let path = tmp.path().join("snap_markers.svg");
        let plot = RegionPlot {
            chrom: "chr1".to_string(),
            start: 10000,
            end: 10500,
            reference: Some(vec![b'A'; 500]),
            genes: Vec::new(),
            samples: Vec::new(),
            markers: vec![
                SnappingMarker {
                    pos: 10020,
                    label: "C>T".to_string(),
                    marker_type: MarkerType::Variant,
                    end_pos: None,
                },
                SnappingMarker {
                    pos: 10100,
                    label: "DEL_1".to_string(),
                    marker_type: MarkerType::StructuralVariant,
                    end_pos: Some(10250),
                },
                SnappingMarker {
                    pos: 10050,
                    label: "ROI_A".to_string(),
                    marker_type: MarkerType::RegionOfInterest,
                    end_pos: Some(10080),
                },
            ],
        };
        render_to_path(&plot, &PlotOptions::default(), &path).unwrap();
        let svg = std::fs::read_to_string(path).unwrap();
        assert!(svg.contains("<svg"));
        assert!(svg.contains("C&gt;T") || svg.contains("C>T"));
        assert!(svg.contains("DEL_1"));
    }
}
fn chrom_match(c1: &str, c2: &str) -> bool {
    normalize_chrom(c1) == normalize_chrom(c2)
}

fn normalize_chrom(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

pub(crate) fn chrom_aliases(chrom: &str) -> Vec<String> {
    let trimmed = chrom.trim();
    if trimmed.is_empty() {
        return Vec::new();
    }

    let lower = trimmed.to_ascii_lowercase();
    let no_chr = lower.strip_prefix("chr").unwrap_or(&lower);
    let canonical_with_chr = format!("chr{}", no_chr);
    let canonical_without_chr = no_chr.to_string();
    let upper_tail = lower
        .strip_prefix("chr")
        .unwrap_or(&lower)
        .to_ascii_uppercase();

    let mut aliases = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for alias in [
        trimmed.to_string(),
        trimmed.to_ascii_uppercase(),
        canonical_with_chr,
        canonical_without_chr,
        lower.clone(),
        upper_tail,
    ] {
        if seen.insert(alias.clone()) {
            aliases.push(alias);
        }
    }

    aliases
}
