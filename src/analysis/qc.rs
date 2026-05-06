use crate::io::bam::match_span;
use noodles::bam;
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
    pub gene_strand: Option<char>,
    pub gene_idx: Option<usize>,
    pub aligned_len: u64,
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
    pub stranded_forward_count: usize,
    pub stranded_reverse_count: usize,
    pub pending_evictions: usize,
    pub inner_distances: Vec<i64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadEndType {
    FivePrime,
    ThreePrime,
}

impl RnaSeqQcSummary {
    pub fn observe_inner_distance(&mut self, dist: i32) {
        self.inner_distances.push(dist as i64);
    }
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
        let stranded_total = self.stranded_forward_count + self.stranded_reverse_count;
        if stranded_total < 100 {
            // Not enough annotated overlap to infer strandedness reliably
            let fr = self.fr_fraction();
            let rf = self.rf_fraction();
            if self.informative_pairs == 0 {
                "undetermined"
            } else if fr >= 0.8 {
                "FR (geom)"
            } else if rf >= 0.8 {
                "RF (geom)"
            } else {
                "mixed (geom)"
            }
        } else {
            let f_frac = self.stranded_forward_count as f64 / stranded_total as f64;
            let r_frac = self.stranded_reverse_count as f64 / stranded_total as f64;
            if f_frac >= 0.8 {
                "FR (stranded-forward)"
            } else if r_frac >= 0.8 {
                "RF (stranded-reverse)"
            } else {
                "unstranded/mixed"
            }
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
        out.push_str(&format!(
            "geometric_fr_fraction\t{:.4}\n",
            self.fr_fraction()
        ));
        out.push_str(&format!(
            "geometric_rf_fraction\t{:.4}\n",
            self.rf_fraction()
        ));
        out.push_str(&format!(
            "other_orientation_fraction\t{:.4}\n",
            self.other_fraction()
        ));
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
        out.push_str(&format!(
            "stranded_forward_count\t{}\n",
            self.stranded_forward_count
        ));
        out.push_str(&format!(
            "stranded_reverse_count\t{}\n",
            self.stranded_reverse_count
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
    pub qc_sample_size: usize,
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
        let _pending_limit = target_pairs.max(2);
        Self {
            qc_sample_size: target_pairs,
            summary: RnaSeqQcSummary {
                requested_pairs: target_pairs,
                ..RnaSeqQcSummary::default()
            },
            pending: HashMap::new(),
            pending_order: VecDeque::new(),
        }
    }

    pub fn needs_more(&self) -> bool {
        self.summary.informative_pairs < self.qc_sample_size
    }

    fn evict_oldest_pending(&mut self) {
        let pending_limit = self.qc_sample_size.max(2);
        while self.pending.len() >= pending_limit {
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

    pub fn observe_read_end(
        &mut self,
        record: &bam::Record,
        chrom: &str,
        dense: &crate::analysis::index::DenseMap,
        _genes: &[crate::models::Gene],
        _end_type: ReadEndType,
    ) {
        if !self.needs_more() {
            return;
        }

        let flags = record.flags();
        let Some(qname) = record.name().map(|n| n.as_ref().to_vec()) else {
            return;
        };
        let Some((start, end)) = match_span(record) else {
            return;
        };

        // Aligned length excluding introns (N)
        let aligned_len = record
            .cigar()
            .iter()
            .map(|result| {
                let op = result.expect("Invalid CIGAR op");
                use noodles::sam::alignment::record::cigar::op::Kind;
                match op.kind() {
                    Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        op.len() as u64
                    }
                    _ => 0,
                }
            })
            .sum();

        // Determine if this read end overlaps a gene on a specific strand
        let mut gene_strand = None;
        let mut gene_idx = None;
        let mid = start + (end - start) / 2;
        use crate::analysis::index::Hits;
        match dense.get_hits(mid) {
            Hits::Single(g_idx) => {
                gene_idx = Some(g_idx as usize);
                gene_strand = Some(_genes[g_idx as usize].representative.strand);
            }
            Hits::Multi(indices) => {
                // If all genes have the same strand, we can use it
                let strands: Vec<char> = indices
                    .iter()
                    .map(|&idx| _genes[idx as usize].representative.strand)
                    .collect();
                if strands.iter().all(|&s| s == strands[0]) {
                    gene_strand = Some(strands[0]);
                    // If multiple genes, just pick the first for index
                    gene_idx = Some(indices[0] as usize);
                }
            }
            Hits::None => {}
        }

        let obs = ReadEndObservation {
            chrom: chrom.to_string(),
            first_match: start,
            last_match: end,
            is_first: flags.is_first_segment(),
            is_reverse: flags.is_reverse_complemented(),
            gene_strand,
            gene_idx,
            aligned_len,
        };

        if let Some((r1, r2)) = self.observe_end(qname, obs) {
            if r1.chrom == r2.chrom {
                // Only use pairs where both hit the same gene for inner distance
                if let (Some(g1), Some(g2)) = (r1.gene_idx, r2.gene_idx) {
                    if g1 == g2 {
                        let gene = &_genes[g1];
                        let p1 = r1.first_match.min(r2.first_match);
                        let p2 = r1.last_match.max(r2.last_match).saturating_sub(1);

                        if let (Some(s1), Some(s2)) = (
                            gene.bin_map.get_spliced_5p(p1),
                            gene.bin_map.get_spliced_5p(p2),
                        ) {
                            let spliced_frag_len = (s1 as i64 - s2 as i64).abs() + 1;
                            let inner =
                                spliced_frag_len - (r1.aligned_len as i64 + r2.aligned_len as i64);

                            // Filter for reasonable inner distance (-1000 to 2000 bp)
                            if inner >= -1000 && inner <= 2000 {
                                self.record_pair(inner, &r1, &r2);
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn record_pair(
        &mut self,
        inner_distance: i64,
        r1: &ReadEndObservation,
        r2: &ReadEndObservation,
    ) {
        if !self.needs_more() {
            return;
        }

        self.summary.inner_distances.push(inner_distance);
        self.summary.informative_pairs += 1;

        // Geometric orientation
        match (r1.is_reverse, r2.is_reverse) {
            (false, true) => self.summary.fr_count += 1,
            (true, false) => self.summary.rf_count += 1,
            _ => self.summary.other_orientation_count += 1,
        }

        // Annotation-relative strandness
        // R1: match(rev, strand=='+') -> RF, else FR
        // R2: match(rev, strand=='+') -> FR, else RF
        if let Some(strand) = r1.gene_strand {
            let matches = r1.is_reverse == (strand == '+');
            if matches {
                self.summary.stranded_reverse_count += 1;
            } else {
                self.summary.stranded_forward_count += 1;
            }
        } else if let Some(strand) = r2.gene_strand {
            let matches = r2.is_reverse == (strand == '+');
            if matches {
                self.summary.stranded_forward_count += 1;
            } else {
                self.summary.stranded_reverse_count += 1;
            }
        }
    }

    pub fn merge_from(&mut self, other: Self) {
        self.summary.aligned_qc_reads += other.summary.aligned_qc_reads;
        self.summary.mtdna_reads += other.summary.mtdna_reads;
        self.summary.rdna_reads += other.summary.rdna_reads;
        // requested_pairs is the target, don't sum it.
        self.summary.informative_pairs += other.summary.informative_pairs;
        self.summary.fr_count += other.summary.fr_count;
        self.summary.rf_count += other.summary.rf_count;
        self.summary.other_orientation_count += other.summary.other_orientation_count;
        self.summary.stranded_forward_count += other.summary.stranded_forward_count;
        self.summary.stranded_reverse_count += other.summary.stranded_reverse_count;
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
            gene_strand: Some('+'),
            gene_idx: Some(0),
            aligned_len: 10,
        };
        let read2 = ReadEndObservation {
            chrom: "chr1".to_string(),
            first_match: 30,
            last_match: 39,
            is_first: false,
            is_reverse: true,
            gene_strand: Some('+'),
            gene_idx: Some(0),
            aligned_len: 10,
        };

        assert!(qc.observe_end(b"r1".to_vec(), read1.clone()).is_none());
        let pair = qc.observe_end(b"r1".to_vec(), read2.clone()).unwrap();
        qc.record_pair(10, &pair.0, &pair.1);

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
            gene_strand: None,
            gene_idx: None,
            aligned_len: 10,
        };

        let _ = qc.observe_end(b"r1".to_vec(), read("chr1", 10));
        let _ = qc.observe_end(b"r2".to_vec(), read("chr1", 20));
        let _ = qc.observe_end(b"r3".to_vec(), read("chr1", 30));

        assert_eq!(qc.summary.pending_evictions, 1);
    }
}
