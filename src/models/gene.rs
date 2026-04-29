use crate::models::Transcript;
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicI32, Ordering};
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

    // Shared Atomic Accumlators (Not serialized)
    #[serde(skip, default = "empty_arc_vec")]
    pub diff_3p: Arc<Vec<AtomicI32>>,
    #[serde(skip, default = "empty_arc_vec")]
    pub diff_percentile: Arc<Vec<AtomicI32>>,
}

fn empty_arc_vec<T>() -> Arc<Vec<T>> {
    Arc::new(Vec::new())
}

impl Gene {
    pub fn new(id: String, name: Option<String>, tx: Transcript, max_3p_dist: usize) -> Self {
        let total_len = tx.total_length as u32;
        let biotype = tx.biotype.clone();
        let bin_map = GeneBinMap::new(&tx);
        let chrom = tx.chrom.clone();

        Self {
            id,
            name,
            chrom,
            biotype,
            total_len,
            representative: tx,
            bin_map: Arc::new(bin_map),
            diff_3p: Arc::new((0..=max_3p_dist).map(|_| AtomicI32::new(0)).collect()),
            diff_percentile: Arc::new((0..101).map(|_| AtomicI32::new(0)).collect()),
        }
    }

    // Removed unused add_coverage_at_genomic_range as it's replaced by thread-local logic.

    /// Apply thread-local accumulated differences to the global atomics.
    pub fn apply_local_updates(&self, d_3p: &[i32], d_pct: &[i32]) {
        for (i, &val) in d_3p.iter().enumerate() {
            if val != 0 && i < self.diff_3p.len() {
                self.diff_3p[i].fetch_add(val, Ordering::Relaxed);
            }
        }
        for (i, &val) in d_pct.iter().enumerate() {
            if val != 0 && i < self.diff_percentile.len() {
                self.diff_percentile[i].fetch_add(val, Ordering::Relaxed);
            }
        }
    }

    pub fn finalize_coverage(&self, diff: &[AtomicI32], out: &mut [u32]) {
        let mut curr = 0i32;
        for i in 0..out.len() {
            curr += diff[i].load(Ordering::Relaxed);
            out[i] = curr.max(0) as u32;
        }
    }

    pub fn reset_coverage(&self) {
        for v in self.diff_3p.iter() {
            v.store(0, Ordering::Relaxed);
        }
        for v in self.diff_percentile.iter() {
            v.store(0, Ordering::Relaxed);
        }
    }

    pub fn reinitialize_atomics(&mut self, max_3p_dist: usize) {
        self.diff_3p = Arc::new((0..=max_3p_dist).map(|_| AtomicI32::new(0)).collect());
        self.diff_percentile = Arc::new((0..101).map(|_| AtomicI32::new(0)).collect());
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
            diff_3p: Arc::clone(&self.diff_3p),
            diff_percentile: Arc::clone(&self.diff_percentile),
        }
    }
}
