use crate::model::ReadModel;

#[derive(Clone, Debug, PartialEq)]
pub struct PlacedRead<'a> {
    pub read: &'a ReadModel,
    pub lane: usize,
}

pub fn place_reads(reads: &[ReadModel]) -> (Vec<PlacedRead<'_>>, usize) {
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|&idx| (reads[idx].start, reads[idx].end, idx));

    let mut lane_ends: Vec<i64> = Vec::new();
    let mut placed = Vec::with_capacity(reads.len());

    for idx in order {
        let read = &reads[idx];
        let lane = lane_ends
            .iter()
            .position(|&end| end <= read.start)
            .unwrap_or_else(|| {
                lane_ends.push(i64::MIN);
                lane_ends.len() - 1
            });

        lane_ends[lane] = read.end.max(read.start + 1);
        placed.push(PlacedRead { read, lane });
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
        let (_, lanes) = place_reads(&reads);
        assert_eq!(lanes, 1);
    }

    #[test]
    fn overlapping_reads_use_extra_lanes() {
        let reads = vec![read("a", 10, 30), read("b", 20, 40)];
        let (_, lanes) = place_reads(&reads);
        assert_eq!(lanes, 2);
    }
}
