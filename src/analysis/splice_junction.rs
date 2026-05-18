use std::collections::{HashSet, HashMap};
use serde::{Serialize, Deserialize};
use crate::io::annotation::ParsedTranscript;

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct JunctionKey {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JunctionNovelty {
    Known,
    PartialNovel,
    Novel,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObservedJunction {
    pub count: u64,
    pub min_hash: u8,
}

#[derive(Clone, Debug)]
pub struct ReferenceJunctions {
    pub known_introns: HashSet<(String, u64, u64)>,
    pub known_starts: HashSet<(String, u64)>,
    pub known_ends: HashSet<(String, u64)>,
}

impl ReferenceJunctions {
    pub fn from_transcripts(transcripts: &[ParsedTranscript]) -> Self {
        let mut known_introns = HashSet::new();
        let mut known_starts = HashSet::new();
        let mut known_ends = HashSet::new();

        for tx in transcripts {
            if tx.exons.len() >= 2 {
                for i in 0..tx.exons.len() - 1 {
                    let start = tx.exons[i].end;
                    let end = tx.exons[i+1].start;
                    if start < end {
                        known_introns.insert((tx.chrom.clone(), start, end));
                        known_starts.insert((tx.chrom.clone(), start));
                        known_ends.insert((tx.chrom.clone(), end));
                    }
                }
            }
        }

        Self {
            known_introns,
            known_starts,
            known_ends,
        }
    }

    pub fn classify(&self, chrom: &str, start: u64, end: u64) -> JunctionNovelty {
        let chrom_str = chrom.to_string();
        if self.known_introns.contains(&(chrom_str.clone(), start, end)) {
            JunctionNovelty::Known
        } else if self.known_starts.contains(&(chrom_str.clone(), start)) || self.known_ends.contains(&(chrom_str, end)) {
            JunctionNovelty::PartialNovel
        } else {
            JunctionNovelty::Novel
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetailedJunction {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub novelty: JunctionNovelty,
    pub read_support: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SpliceJunctionMetrics {
    pub unique_known: usize,
    pub unique_partial_novel: usize,
    pub unique_novel: usize,
    pub total_known_reads: u64,
    pub total_partial_novel_reads: u64,
    pub total_novel_reads: u64,
    pub saturation_known: Vec<usize>,
    pub saturation_partial_novel: Vec<usize>,
    pub saturation_novel: Vec<usize>,
    pub details: Vec<DetailedJunction>,
}

pub fn compute_metrics(
    junctions: &HashMap<JunctionKey, ObservedJunction>,
    ref_juncs: &ReferenceJunctions,
) -> SpliceJunctionMetrics {
    let mut unique_known = 0;
    let mut unique_partial_novel = 0;
    let mut unique_novel = 0;
    let mut total_known_reads = 0;
    let mut total_partial_novel_reads = 0;
    let mut total_novel_reads = 0;

    let mut sat_known = vec![0; 10];
    let mut sat_partial = vec![0; 10];
    let mut sat_novel = vec![0; 10];
    let mut details = Vec::new();

    for (key, obs) in junctions {
        let novelty = ref_juncs.classify(&key.chrom, key.start, key.end);
        match novelty {
            JunctionNovelty::Known => {
                unique_known += 1;
                total_known_reads += obs.count;
            }
            JunctionNovelty::PartialNovel => {
                unique_partial_novel += 1;
                total_partial_novel_reads += obs.count;
            }
            JunctionNovelty::Novel => {
                unique_novel += 1;
                total_novel_reads += obs.count;
            }
        }

        for bin in 0..10 {
            let threshold = (bin + 1) * 10;
            if obs.min_hash < threshold as u8 {
                match novelty {
                    JunctionNovelty::Known => sat_known[bin] += 1,
                    JunctionNovelty::PartialNovel => sat_partial[bin] += 1,
                    JunctionNovelty::Novel => sat_novel[bin] += 1,
                }
            }
        }

        details.push(DetailedJunction {
            chrom: key.chrom.clone(),
            start: key.start,
            end: key.end,
            strand: key.strand,
            novelty,
            read_support: obs.count,
        });
    }

    // Sort details genomic position for deterministic clean output
    details.sort_by(|a, b| {
        a.chrom.cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    SpliceJunctionMetrics {
        unique_known,
        unique_partial_novel,
        unique_novel,
        total_known_reads,
        total_partial_novel_reads,
        total_novel_reads,
        saturation_known: sat_known,
        saturation_partial_novel: sat_partial,
        saturation_novel: sat_novel,
        details,
    }
}
