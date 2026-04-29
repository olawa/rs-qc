use crate::analysis::read_distribution::RegionType;
use crate::io::annotation::ParsedTranscript;
use crate::models::Exon;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Arc;

const REGION_ORDER: [RegionType; 11] = [
    RegionType::CdsExon,
    RegionType::Utr5Exon,
    RegionType::Utr3Exon,
    RegionType::Exon,
    RegionType::Intron,
    RegionType::TssUp1kb,
    RegionType::TssUp5kb,
    RegionType::TssUp10kb,
    RegionType::TesDown1kb,
    RegionType::TesDown5kb,
    RegionType::TesDown10kb,
];

#[derive(Debug, Clone, Copy)]
struct FeatureEvent {
    pos: u64,
    region: RegionType,
    delta: i32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureInterval {
    pub start: u64,
    pub end: u64,
    pub region: RegionType,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromFeatureIndex {
    pub intervals: Arc<[FeatureInterval]>,
}

#[derive(Debug, Default, Clone)]
pub struct FeatureCursor {
    idx: usize,
}

#[derive(Debug, Clone)]
pub struct FeatureIndex {
    pub chroms: HashMap<String, Arc<ChromFeatureIndex>>,
    pub feature_sizes: HashMap<RegionType, u64>,
}

impl FeatureIndex {
    pub fn build(transcripts: &[ParsedTranscript]) -> Self {
        let mut chrom_events: HashMap<String, Vec<FeatureEvent>> = HashMap::new();

        for tx in transcripts {
            let events = chrom_events.entry(tx.chrom.clone()).or_default();
            emit_transcript_events(tx, events);
        }

        let mut chroms = HashMap::new();
        let mut feature_sizes: HashMap<RegionType, u64> = HashMap::new();

        for (chrom, mut events) in chrom_events {
            if events.is_empty() {
                continue;
            }

            events.sort_by(|a, b| a.pos.cmp(&b.pos));

            let mut active = [0i32; 12];
            let mut intervals = Vec::<FeatureInterval>::new();
            let mut prev_pos: Option<u64> = None;
            let mut i = 0usize;

            while i < events.len() {
                let pos = events[i].pos;

                if let Some(prev) = prev_pos {
                    if prev < pos {
                        let region = pick_region(&active);
                        if region != RegionType::Intergenic {
                            let len = pos - prev;
                            *feature_sizes.entry(region).or_insert(0) += len;
                            push_interval(&mut intervals, prev, pos, region);
                        }
                    }
                }

                while i < events.len() && events[i].pos == pos {
                    let region_idx = events[i].region as usize;
                    active[region_idx] += events[i].delta;
                    i += 1;
                }
                prev_pos = Some(pos);
            }

            if !intervals.is_empty() {
                chroms.insert(
                    chrom,
                    Arc::new(ChromFeatureIndex {
                        intervals: intervals.into(),
                    }),
                );
            }
        }

        Self {
            chroms,
            feature_sizes,
        }
    }
}

impl ChromFeatureIndex {
    pub fn cursor_at(&self, pos: u64) -> FeatureCursor {
        FeatureCursor {
            idx: self.intervals.partition_point(|iv| iv.end <= pos),
        }
    }

    pub fn classify(&self, pos: u64, cursor: &mut FeatureCursor) -> RegionType {
        while cursor.idx < self.intervals.len() && self.intervals[cursor.idx].end <= pos {
            cursor.idx += 1;
        }

        if let Some(interval) = self.intervals.get(cursor.idx) {
            if interval.start <= pos && pos < interval.end {
                return interval.region;
            }
        }

        RegionType::Intergenic
    }
}

fn emit_transcript_events(tx: &ParsedTranscript, events: &mut Vec<FeatureEvent>) {
    let mut exons = tx.exons.clone();
    exons.sort_by_key(|e| e.start);

    for exon in &exons {
        emit_exon_events(tx, exon, events);
    }

    for i in 0..exons.len().saturating_sub(1) {
        let intron_s = exons[i].end;
        let intron_e = exons[i + 1].start;
        if intron_s < intron_e {
            push_event(events, intron_s, intron_e, RegionType::Intron);
        }
    }

    if let (Some(tss), Some(tes)) = (tx_tss(tx, &exons), tx_tes(tx, &exons)) {
        emit_terminal_windows(tx.strand, tss, tes, events);
    }
}

fn emit_exon_events(tx: &ParsedTranscript, exon: &Exon, events: &mut Vec<FeatureEvent>) {
    if let (Some(cds_s), Some(cds_e)) = (tx.cds_start, tx.cds_end) {
        let overlap_s = exon.start.max(cds_s);
        let overlap_e = exon.end.min(cds_e);
        if overlap_s < overlap_e {
            push_event(events, overlap_s, overlap_e, RegionType::CdsExon);
        }

        if tx.strand == '+' {
            let utr5_e = exon.end.min(cds_s);
            if exon.start < utr5_e {
                push_event(events, exon.start, utr5_e, RegionType::Utr5Exon);
            }

            let utr3_s = exon.start.max(cds_e);
            if utr3_s < exon.end {
                push_event(events, utr3_s, exon.end, RegionType::Utr3Exon);
            }
        } else {
            let utr5_s = exon.start.max(cds_e);
            if utr5_s < exon.end {
                push_event(events, utr5_s, exon.end, RegionType::Utr5Exon);
            }

            let utr3_e = exon.end.min(cds_s);
            if exon.start < utr3_e {
                push_event(events, exon.start, utr3_e, RegionType::Utr3Exon);
            }
        }
    } else {
        push_event(events, exon.start, exon.end, RegionType::Exon);
    }
}

fn emit_terminal_windows(strand: char, tss: u64, tes: u64, events: &mut Vec<FeatureEvent>) {
    if strand == '+' {
        push_event(
            events,
            tss.saturating_sub(10_000),
            tss.saturating_sub(5_000),
            RegionType::TssUp10kb,
        );
        push_event(
            events,
            tss.saturating_sub(5_000),
            tss.saturating_sub(1_000),
            RegionType::TssUp5kb,
        );
        push_event(events, tss.saturating_sub(1_000), tss, RegionType::TssUp1kb);

        push_event(
            events,
            tes,
            tes.saturating_add(1_000),
            RegionType::TesDown1kb,
        );
        push_event(
            events,
            tes.saturating_add(1_000),
            tes.saturating_add(5_000),
            RegionType::TesDown5kb,
        );
        push_event(
            events,
            tes.saturating_add(5_000),
            tes.saturating_add(10_000),
            RegionType::TesDown10kb,
        );
    } else {
        push_event(
            events,
            tss.saturating_add(5_000),
            tss.saturating_add(10_000),
            RegionType::TssUp10kb,
        );
        push_event(
            events,
            tss.saturating_add(1_000),
            tss.saturating_add(5_000),
            RegionType::TssUp5kb,
        );
        push_event(events, tss, tss.saturating_add(1_000), RegionType::TssUp1kb);

        push_event(
            events,
            tes.saturating_sub(1_000),
            tes,
            RegionType::TesDown1kb,
        );
        push_event(
            events,
            tes.saturating_sub(5_000),
            tes.saturating_sub(1_000),
            RegionType::TesDown5kb,
        );
        push_event(
            events,
            tes.saturating_sub(10_000),
            tes.saturating_sub(5_000),
            RegionType::TesDown10kb,
        );
    }
}

fn tx_tss(tx: &ParsedTranscript, exons: &[Exon]) -> Option<u64> {
    if exons.is_empty() {
        return None;
    }
    Some(if tx.strand == '+' {
        exons.first()?.start
    } else {
        exons.last()?.end
    })
}

fn tx_tes(tx: &ParsedTranscript, exons: &[Exon]) -> Option<u64> {
    if exons.is_empty() {
        return None;
    }
    Some(if tx.strand == '+' {
        exons.last()?.end
    } else {
        exons.first()?.start
    })
}

fn push_event(events: &mut Vec<FeatureEvent>, start: u64, end: u64, region: RegionType) {
    if start >= end {
        return;
    }
    events.push(FeatureEvent {
        pos: start,
        region,
        delta: 1,
    });
    events.push(FeatureEvent {
        pos: end,
        region,
        delta: -1,
    });
}

fn pick_region(active: &[i32; 12]) -> RegionType {
    for region in REGION_ORDER {
        if active[region as usize] > 0 {
            return region;
        }
    }
    RegionType::Intergenic
}

fn push_interval(intervals: &mut Vec<FeatureInterval>, start: u64, end: u64, region: RegionType) {
    if start >= end {
        return;
    }
    if let Some(last) = intervals.last_mut() {
        if last.region == region && last.end == start {
            last.end = end;
            return;
        }
    }

    intervals.push(FeatureInterval { start, end, region });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::annotation::IdResolutionStrategy;

    fn transcript(
        gene_id: &str,
        chrom: &str,
        strand: char,
        exons: Vec<(u64, u64)>,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
    ) -> ParsedTranscript {
        ParsedTranscript {
            gene_id: gene_id.to_string(),
            transcript_id: "tx".to_string(),
            gene_name: None,
            biotype: None,
            chrom: chrom.to_string(),
            strand,
            exons: exons
                .into_iter()
                .map(|(start, end)| Exon { start, end })
                .collect(),
            cds_start,
            cds_end,
            strategy: IdResolutionStrategy::GtfExplicit,
        }
    }

    #[test]
    fn classifies_cds_utr_and_intron() {
        let tx = transcript(
            "g1",
            "chr1",
            '+',
            vec![(100, 200), (300, 400)],
            Some(120),
            Some(350),
        );
        let idx = FeatureIndex::build(&[tx]);
        let chrom = idx.chroms.get("chr1").unwrap();
        let mut cursor = chrom.cursor_at(0);

        assert_eq!(chrom.classify(110, &mut cursor), RegionType::Utr5Exon);
        assert_eq!(chrom.classify(130, &mut cursor), RegionType::CdsExon);
        assert_eq!(chrom.classify(210, &mut cursor), RegionType::Intron);
        assert_eq!(chrom.classify(360, &mut cursor), RegionType::Utr3Exon);
    }
}
