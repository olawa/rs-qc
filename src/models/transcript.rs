use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub struct Exon {
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Transcript {
    pub id: String,
    pub chrom: String,
    pub strand: char,
    pub biotype: Option<String>,
    pub exons: Vec<Exon>, // Sorted by genomic position
    pub total_length: u64,
    pub cds_start: Option<u64>, // Genomic 0-based
    pub cds_end: Option<u64>,   // Genomic 0-based, half-open
}

impl Transcript {
    pub fn new(
        id: String,
        chrom: String,
        strand: char,
        biotype: Option<String>,
        mut exons: Vec<Exon>,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
    ) -> Self {
        exons.sort_by_key(|e| e.start);
        let total_length = exons.iter().map(|e| e.end - e.start).sum();
        Self {
            id,
            chrom,
            strand,
            biotype,
            exons,
            total_length,
            cds_start,
            cds_end,
        }
    }

    /// Returns the 5' terminal genomic coordinate.
    pub fn five_prime_end(&self) -> u64 {
        if self.strand == '+' {
            self.exons.first().unwrap().start
        } else {
            self.exons.last().unwrap().end
        }
    }

    /// Returns the 3' terminal genomic coordinate.
    pub fn three_prime_end(&self) -> u64 {
        if self.strand == '+' {
            self.exons.last().unwrap().end
        } else {
            self.exons.first().unwrap().start
        }
    }

    /// Returns (min_start, max_end) of all exons.
    pub fn genomic_span(&self) -> (u64, u64) {
        if self.exons.is_empty() {
            return (0, 0);
        }
        (
            self.exons.first().unwrap().start,
            self.exons.last().unwrap().end,
        )
    }

    pub fn genomic_range_to_spliced_3p(&self, start: u64, end: u64) -> Vec<(u64, u64)> {
        let mut results = Vec::new();
        // Simply reuse the 5' logic and flip it
        let intervals_5p = self.genomic_range_to_spliced_5p(start, end);
        for (tx_start, tx_end) in intervals_5p {
            // Biological start 0 is 5'.
            // Distance from 3' = total_length - 1 - dist_from_5p
            // Range [s, e) from 5' maps to [L-e, L-s) from 3'
            let s_3p = self.total_length.saturating_sub(tx_end);
            let e_3p = self.total_length.saturating_sub(tx_start);
            results.push((s_3p, e_3p));
        }
        results
    }

    pub fn genomic_range_to_spliced_5p(&self, start: u64, end: u64) -> Vec<(u64, u64)> {
        let mut results = Vec::new();
        if self.strand == '+' {
            let mut offset = 0;
            for exon in &self.exons {
                let overlap_start = start.max(exon.start);
                let overlap_end = end.min(exon.end);
                if overlap_start < overlap_end {
                    results.push((
                        offset + (overlap_start - exon.start),
                        offset + (overlap_end - exon.start),
                    ));
                }
                offset += exon.end - exon.start;
            }
        } else {
            let mut offset = 0;
            for exon in self.exons.iter().rev() {
                let overlap_start = start.max(exon.start);
                let overlap_end = end.min(exon.end);
                if overlap_start < overlap_end {
                    // For - strand, genomic start of exon is biological end
                    results.push((
                        offset + (exon.end - overlap_end),
                        offset + (exon.end - overlap_start),
                    ));
                }
                offset += exon.end - exon.start;
            }
        }
        results
    }

    /// Maps a 5' spliced position (0..total_length-1) to a genomic coordinate.
    pub fn spliced_to_genomic(&self, spliced_pos: u64) -> Option<u64> {
        if spliced_pos >= self.total_length {
            return None;
        }

        if self.strand == '+' {
            let mut offset = 0;
            for exon in &self.exons {
                let len = exon.end - exon.start;
                if spliced_pos >= offset && spliced_pos < offset + len {
                    return Some(exon.start + (spliced_pos - offset));
                }
                offset += len;
            }
        } else {
            let mut offset = 0;
            for exon in self.exons.iter().rev() {
                let len = exon.end - exon.start;
                if spliced_pos >= offset && spliced_pos < offset + len {
                    // For - strand, biological start of transcript (5') is genomic END of last exon
                    return Some(exon.end - 1 - (spliced_pos - offset));
                }
                offset += len;
            }
        }
        None
    }
}
