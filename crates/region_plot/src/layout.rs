use crate::model::{GeneModel, ReadModel};

#[derive(Clone, Debug, PartialEq)]
pub struct PlacedRead<'a> {
    pub read: &'a ReadModel,
    pub lane: usize,
}

pub fn place_reads(reads: &[ReadModel], squash: bool) -> (Vec<PlacedRead<'_>>, usize) {
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|&idx| (reads[idx].start, reads[idx].end, idx));

    let mut lane_ends: Vec<i64> = Vec::new();
    let mut placed = Vec::with_capacity(reads.len());

    let gap = 0;

    for idx in order {
        let read = &reads[idx];
        let lane = lane_ends
            .iter()
            .position(|&end| end + gap <= read.start)
            .unwrap_or_else(|| {
                lane_ends.push(i64::MIN);
                lane_ends.len() - 1
            });

        lane_ends[lane] = read.end.max(read.start + 1);
        placed.push(PlacedRead { read, lane });
    }

    let mut total_lanes = lane_ends.len();
    if squash {
        // Limit lanes to something reasonable if we have too many.
        total_lanes = total_lanes.min(500);
        placed.retain(|p| p.lane < 500);
    }
    (placed, total_lanes)
}

pub struct PlacedGene<'a> {
    pub gene: &'a GeneModel,
    pub lane: usize,
}

pub fn place_genes(genes: &[GeneModel]) -> (Vec<PlacedGene<'_>>, usize) {
    let mut order: Vec<usize> = (0..genes.len()).collect();
    order.sort_by_key(|&idx| (genes[idx].start, genes[idx].end, idx));

    let mut lane_ends: Vec<i64> = Vec::new();
    let mut placed = Vec::with_capacity(genes.len());

    // GTF genes usually need more horizontal padding between names
    let gap = 500;

    for idx in order {
        let gene = &genes[idx];
        let lane = lane_ends
            .iter()
            .position(|&end| end + gap <= gene.start)
            .unwrap_or_else(|| {
                lane_ends.push(i64::MIN);
                lane_ends.len() - 1
            });

        lane_ends[lane] = gene.end;
        placed.push(PlacedGene { gene, lane });
    }

    let total_lanes = lane_ends.len();
    (placed, total_lanes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::ReadSegment;

    fn read(name: &str, start: i64, end: i64) -> ReadModel {
        ReadModel {
            name: name.to_string(),
            start,
            end,
            is_reverse: false,
            mapq: 60,
            segments: vec![ReadSegment::Match {
                ref_start: start,
                len: end - start,
                query_start: 0,
            }],
            bases: None,
            qualities: None,
            haplotype: None,
            modifications: Vec::new(),
        }
    }

    #[test]
    fn non_overlapping_reads_share_lanes() {
        let reads = vec![read("a", 10, 20), read("b", 20, 30)];
        let (_, lanes) = place_reads(&reads, false);
        assert_eq!(lanes, 1);
    }

    #[test]
    fn overlapping_reads_use_extra_lanes() {
        let reads = vec![read("a", 10, 30), read("b", 20, 40)];
        let (_, lanes) = place_reads(&reads, false);
        assert_eq!(lanes, 2);
    }

    #[test]
    fn squash_mode_still_separates_overlaps() {
        let reads = vec![read("a", 10, 30), read("b", 20, 40)];
        let (_, lanes) = place_reads(&reads, true);
        assert_eq!(lanes, 2);
    }
}
