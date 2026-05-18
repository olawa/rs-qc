use crate::models::Transcript;
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Arc;

/// A pre-computed mapping from genomic offset to biological 5' spliced position.
/// This enables O(1) bin lookups in the hot loop.
#[derive(Debug, Serialize, Deserialize)]
pub struct GeneBinMap {
    pub offset: u64,
    pub bins: Vec<u32>, // Offset (pos - start) -> 0-based spliced pos (from 5')
}

impl GeneBinMap {
    pub fn new(tx: &Transcript) -> Self {
        let (span_start, span_end) = tx.genomic_span();
        let span = (span_end - span_start) as usize;

        // Use u32::MAX as sentinel for non-exonic positions
        let mut bins = vec![u32::MAX; span];

        if tx.strand == '+' {
            let mut spliced_pos = 0u32;
            for exon in &tx.exons {
                for pos in exon.start..exon.end {
                    let offset = (pos - span_start) as usize;
                    if offset < span {
                        bins[offset] = spliced_pos;
                        spliced_pos += 1;
                    }
                }
            }
        } else {
            let mut spliced_pos = 0u32;
            for exon in tx.exons.iter().rev() {
                for pos in (exon.start..exon.end).rev() {
                    let offset = (pos - span_start) as usize;
                    if offset < span {
                        bins[offset] = spliced_pos;
                        spliced_pos += 1;
                    }
                }
            }
        }

        Self {
            offset: span_start,
            bins,
        }
    }

    #[inline(always)]
    pub fn get_spliced_5p(&self, genomic_pos: u64) -> Option<u32> {
        if genomic_pos < self.offset {
            return None;
        }
        let l_offset = (genomic_pos - self.offset) as usize;
        if l_offset < self.bins.len() {
            let val = self.bins[l_offset];
            if val != u32::MAX {
                return Some(val);
            }
        }
        None
    }
}

#[derive(Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: Option<String>,
    pub chrom: String,
    pub biotype: Option<String>,
    pub total_len: u32,
    pub representative: Transcript,
    pub bin_map: Arc<GeneBinMap>,

    // Shared Atomic Accumulators (Direct counts, not diffs)
    #[serde(skip, default = "empty_arc_vec")]
    pub counts_3p: Arc<Vec<AtomicU32>>,
    #[serde(skip, default = "empty_arc_vec")]
    pub counts_percentile: Arc<Vec<AtomicU32>>,
    pub num_isoforms: usize,
}

fn empty_arc_vec<T>() -> Arc<Vec<T>> {
    Arc::new(Vec::new())
}

impl Gene {
    pub fn new(
        id: String,
        name: Option<String>,
        tx: Transcript,
        max_3p_dist: usize,
        bin_size: usize,
    ) -> Self {
        let total_len = tx.total_length as u32;
        let biotype = tx.biotype.clone();
        let bin_map = GeneBinMap::new(&tx);
        let chrom = tx.chrom.clone();
        let n_3p_bins = (max_3p_dist + bin_size - 1) / bin_size;

        Self {
            id,
            name,
            chrom,
            biotype,
            total_len,
            representative: tx,
            bin_map: Arc::new(bin_map),
            counts_3p: Arc::new((0..n_3p_bins).map(|_| AtomicU32::new(0)).collect()),
            counts_percentile: Arc::new((0..100).map(|_| AtomicU32::new(0)).collect()),
            num_isoforms: 1,
        }
    }

    pub fn with_num_isoforms(mut self, num_isoforms: usize) -> Self {
        self.num_isoforms = num_isoforms;
        self
    }

    pub fn add_percentile(&self, idx: usize, delta: u32) {
        if idx < self.counts_percentile.len() {
            self.counts_percentile[idx].fetch_add(delta, Ordering::Relaxed);
        }
    }

    pub fn add_3p(&self, idx: usize, delta: u32) {
        if idx < self.counts_3p.len() {
            self.counts_3p[idx].fetch_add(delta, Ordering::Relaxed);
        }
    }

    pub fn reset_coverage(&self) {
        for v in self.counts_3p.iter() {
            v.store(0, Ordering::Relaxed);
        }
        for v in self.counts_percentile.iter() {
            v.store(0, Ordering::Relaxed);
        }
    }

    pub fn reinitialize_atomics(&mut self, max_3p_dist: usize, bin_size: usize) {
        let n_3p_bins = (max_3p_dist + bin_size - 1) / bin_size;
        self.counts_3p = Arc::new((0..n_3p_bins).map(|_| AtomicU32::new(0)).collect());
        self.counts_percentile = Arc::new((0..100).map(|_| AtomicU32::new(0)).collect());
    }
}

impl Clone for Gene {
    fn clone(&self) -> Self {
        Self {
            id: self.id.clone(),
            name: self.name.clone(),
            biotype: self.biotype.clone(),
            chrom: self.chrom.clone(),
            total_len: self.total_len,
            representative: self.representative.clone(),
            bin_map: Arc::clone(&self.bin_map),
            counts_3p: Arc::clone(&self.counts_3p),
            counts_percentile: Arc::clone(&self.counts_percentile),
            num_isoforms: self.num_isoforms,
        }
    }
}
