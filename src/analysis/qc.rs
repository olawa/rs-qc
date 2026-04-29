use serde::Serialize;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs::File;
use std::io::Write;

#[derive(Debug, Clone)]
pub struct ReadEndObservation {
    pub chrom: String,
    pub first_match: u64,
    pub last_match: u64,
    pub is_first: bool,
    pub is_reverse: bool,
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct RnaSeqQcSummary {
    pub aligned_qc_reads: u64,
    pub mtdna_reads: u64,
    pub rdna_reads: u64,
    pub requested_pairs: usize,
    pub informative_pairs: usize,
    pub fr_count: usize,
    pub rf_count: usize,
    pub other_orientation_count: usize,
    pub pending_evictions: usize,
    pub inner_distances: Vec<i64>,
}

impl RnaSeqQcSummary {
    pub fn mtdna_fraction(&self) -> f64 {
        frac_u64(self.mtdna_reads, self.aligned_qc_reads)
    }

    pub fn rdna_fraction(&self) -> f64 {
        frac_u64(self.rdna_reads, self.aligned_qc_reads)
    }

    pub fn fr_fraction(&self) -> f64 {
        frac(self.fr_count, self.informative_pairs)
    }

    pub fn rf_fraction(&self) -> f64 {
        frac(self.rf_count, self.informative_pairs)
    }

    pub fn other_fraction(&self) -> f64 {
        frac(self.other_orientation_count, self.informative_pairs)
    }

    pub fn inferred_strandness(&self) -> &'static str {
        let fr = self.fr_fraction();
        let rf = self.rf_fraction();
        if self.informative_pairs == 0 {
            "undetermined"
        } else if fr >= 0.8 {
            "FR"
        } else if rf >= 0.8 {
            "RF"
        } else {
            "unstranded_or_mixed"
        }
    }

    pub fn inner_distance_mean(&self) -> Option<f64> {
        if self.inner_distances.is_empty() {
            return None;
        }
        let sum: i64 = self.inner_distances.iter().sum();
        Some(sum as f64 / self.inner_distances.len() as f64)
    }

    pub fn inner_distance_median(&self) -> Option<f64> {
        if self.inner_distances.is_empty() {
            return None;
        }
        let mut values = self.inner_distances.clone();
        values.sort_unstable();
        let mid = values.len() / 2;
        if values.len() % 2 == 0 {
            Some((values[mid - 1] as f64 + values[mid] as f64) / 2.0)
        } else {
            Some(values[mid] as f64)
        }
    }

    pub fn inner_distance_histogram(&self) -> BTreeMap<i64, usize> {
        let mut hist = BTreeMap::new();
        for &dist in &self.inner_distances {
            *hist.entry(dist).or_insert(0) += 1;
        }
        hist
    }

    pub fn clipped_inner_distance_series(&self, min_dist: i64, max_dist: i64) -> Vec<f64> {
        let len = (max_dist - min_dist + 1).max(0) as usize;
        let mut counts = vec![0.0; len];
        for &dist in &self.inner_distances {
            if dist < min_dist || dist > max_dist {
                continue;
            }
            counts[(dist - min_dist) as usize] += 1.0;
        }
        counts
    }

    pub fn summary_text(&self) -> String {
        let mut out = String::new();
        out.push_str("metric\tvalue\n");
        out.push_str(&format!("aligned_qc_reads\t{}\n", self.aligned_qc_reads));
        out.push_str(&format!("mtdna_reads\t{}\n", self.mtdna_reads));
        out.push_str(&format!("mtdna_fraction\t{:.4}\n", self.mtdna_fraction()));
        out.push_str(&format!("rdna_reads\t{}\n", self.rdna_reads));
        out.push_str(&format!("rdna_fraction\t{:.4}\n", self.rdna_fraction()));
        out.push_str(&format!("requested_pairs\t{}\n", self.requested_pairs));
        out.push_str(&format!("informative_pairs\t{}\n", self.informative_pairs));
        out.push_str(&format!(
            "inferred_strandness\t{}\n",
            self.inferred_strandness()
        ));
        out.push_str(&format!("fr_fraction\t{:.4}\n", self.fr_fraction()));
        out.push_str(&format!("rf_fraction\t{:.4}\n", self.rf_fraction()));
        out.push_str(&format!("other_fraction\t{:.4}\n", self.other_fraction()));
        out.push_str(&format!("pending_evictions\t{}\n", self.pending_evictions));
        out.push_str(&format!(
            "inner_distance_mean\t{}\n",
            self.inner_distance_mean()
                .map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".to_string())
        ));
        out.push_str(&format!(
            "inner_distance_median\t{}\n",
            self.inner_distance_median()
                .map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".to_string())
        ));
        out
    }

    pub fn write_summary_file(&self, path: &str) -> std::io::Result<()> {
        std::fs::write(path, self.summary_text())
    }

    pub fn write_inner_distance_histogram(&self, path: &str) -> std::io::Result<()> {
        let hist = self.inner_distance_histogram();
        let mut file = File::create(path)?;
        writeln!(file, "distance\tcount")?;
        for (dist, count) in hist {
            writeln!(file, "{}\t{}", dist, count)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct InlineQcState {
    target_pairs: usize,
    pending_limit: usize,
    pub summary: RnaSeqQcSummary,
    pending: HashMap<Vec<u8>, ReadEndObservation>,
    pending_order: VecDeque<Vec<u8>>,
}

impl Default for InlineQcState {
    fn default() -> Self {
        Self::new(100_000)
    }
}

impl InlineQcState {
    pub fn new(target_pairs: usize) -> Self {
        let target_pairs = target_pairs.max(1);
        let pending_limit = target_pairs.max(2);
        Self {
            target_pairs,
            pending_limit,
            summary: RnaSeqQcSummary {
                requested_pairs: target_pairs,
                ..RnaSeqQcSummary::default()
            },
            pending: HashMap::new(),
            pending_order: VecDeque::new(),
        }
    }

    pub fn needs_more(&self) -> bool {
        self.summary.informative_pairs < self.target_pairs
    }

    fn evict_oldest_pending(&mut self) {
        while self.pending.len() >= self.pending_limit {
            let Some(qname) = self.pending_order.pop_front() else {
                break;
            };
            if self.pending.remove(&qname).is_some() {
                self.summary.pending_evictions += 1;
                break;
            }
        }
    }

    pub fn observe_end(
        &mut self,
        qname: Vec<u8>,
        mate: ReadEndObservation,
    ) -> Option<(ReadEndObservation, ReadEndObservation)> {
        if !self.needs_more() {
            return None;
        }

        if let Some(prev) = self.pending.remove(&qname) {
            if prev.is_first == mate.is_first {
                return None;
            }

            if prev.is_first {
                Some((prev, mate))
            } else {
                Some((mate, prev))
            }
        } else {
            self.evict_oldest_pending();
            self.pending_order.push_back(qname.clone());
            self.pending.insert(qname, mate);
            None
        }
    }

    pub fn record_pair(
        &mut self,
        inner_distance: i64,
        read1_is_reverse: bool,
        read2_is_reverse: bool,
    ) {
        if !self.needs_more() {
            return;
        }

        self.summary.inner_distances.push(inner_distance);
        self.summary.informative_pairs += 1;

        match (read1_is_reverse, read2_is_reverse) {
            (false, true) => self.summary.fr_count += 1,
            (true, false) => self.summary.rf_count += 1,
            _ => self.summary.other_orientation_count += 1,
        }
    }

    pub fn merge_from(&mut self, other: Self) {
        self.summary.requested_pairs += other.summary.requested_pairs;
        self.summary.informative_pairs += other.summary.informative_pairs;
        self.summary.fr_count += other.summary.fr_count;
        self.summary.rf_count += other.summary.rf_count;
        self.summary.other_orientation_count += other.summary.other_orientation_count;
        self.summary.pending_evictions += other.summary.pending_evictions;
        self.summary
            .inner_distances
            .extend(other.summary.inner_distances);
    }
}

fn frac(n: usize, d: usize) -> f64 {
    if d == 0 {
        0.0
    } else {
        n as f64 / d as f64
    }
}

fn frac_u64(n: u64, d: u64) -> f64 {
    if d == 0 {
        0.0
    } else {
        n as f64 / d as f64
    }
}

#[cfg(test)]
mod tests {
    use super::{InlineQcState, ReadEndObservation};

    #[test]
    fn pairs_are_classified_and_scored() {
        let mut qc = InlineQcState::new(10);
        let read1 = ReadEndObservation {
            chrom: "chr1".to_string(),
            first_match: 10,
            last_match: 19,
            is_first: true,
            is_reverse: false,
        };
        let read2 = ReadEndObservation {
            chrom: "chr1".to_string(),
            first_match: 30,
            last_match: 39,
            is_first: false,
            is_reverse: true,
        };

        assert!(qc.observe_end(b"r1".to_vec(), read1).is_none());
        let pair = qc.observe_end(b"r1".to_vec(), read2).unwrap();
        qc.record_pair(10, pair.0.is_reverse, pair.1.is_reverse);

        assert_eq!(qc.summary.informative_pairs, 1);
        assert_eq!(qc.summary.fr_count, 1);
        assert_eq!(qc.summary.inner_distances, vec![10]);
    }

    #[test]
    fn oldest_unmatched_reads_are_evicted() {
        let mut qc = InlineQcState::new(2);
        let read = |chrom: &str, pos: u64| ReadEndObservation {
            chrom: chrom.to_string(),
            first_match: pos,
            last_match: pos + 10,
            is_first: true,
            is_reverse: false,
        };

        let _ = qc.observe_end(b"r1".to_vec(), read("chr1", 10));
        let _ = qc.observe_end(b"r2".to_vec(), read("chr1", 20));
        let _ = qc.observe_end(b"r3".to_vec(), read("chr1", 30));

        assert_eq!(qc.summary.pending_evictions, 1);
    }
}
