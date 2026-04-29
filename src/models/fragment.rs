use crate::models::Exon;

pub struct Fragment {
    pub chrom: String,
    pub aligned_blocks: Vec<Exon>, // Genomic segments representing the physical molecule
}

impl Fragment {
    /// Merges multiple aligned blocks into a minimal set of non-overlapping intervals.
    /// Deduplicates overlaps but preserves intronic gaps.
    pub fn new(chrom: String, mut blocks: Vec<Exon>) -> Self {
        if blocks.is_empty() {
            return Self {
                chrom,
                aligned_blocks: Vec::new(),
            };
        }

        blocks.sort_by_key(|e| e.start);

        let mut merged = Vec::with_capacity(blocks.len());
        let mut current = blocks[0].clone();

        for next in blocks.into_iter().skip(1) {
            if next.start < current.end {
                // Overlap found
                if next.end > current.end {
                    current.end = next.end;
                }
            } else {
                // Gap found (could be an intron or just space between mate pairs)
                merged.push(current);
                current = next;
            }
        }
        merged.push(current);

        Self {
            chrom,
            aligned_blocks: merged,
        }
    }
}
