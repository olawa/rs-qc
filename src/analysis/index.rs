use crate::analysis::read_distribution::{build_feature_maps, FeatureMap, RegionType};
use crate::models::{normalize_chrom, Gene};
use std::collections::HashMap;
use std::sync::Arc;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct AnnotationIndex {
    pub version: u32,
    pub genes: Arc<Vec<Gene>>,
    pub feature_sizes: HashMap<RegionType, u64>,
    // We will build dense maps per-chromosome as needed to save memory
    pub chrom_spans: HashMap<String, (u64, u64)>,
    pub genes_by_chrom: HashMap<String, Vec<usize>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DenseMap {
    pub offset: u64,
    pub map: Vec<u32>, // 0 = None, N < 0x80000000 = genes[N-1] + 1, N >= 0x80000000 = multi_map index
    pub multi_data: Vec<u32>,
    pub multi_offsets: Vec<u32>,
}

pub enum Hits<'a> {
    None,
    Single(u32),
    Multi(&'a [u32]),
}

impl AnnotationIndex {
    pub fn new(genes: Vec<Gene>, calculate_features: bool) -> Self {
        let genes = Arc::new(genes);
        let mut chrom_spans: HashMap<String, (u64, u64)> = HashMap::new();
        let mut chrom_max: HashMap<String, u64> = HashMap::new();
        let mut genes_by_chrom: HashMap<String, Vec<usize>> = HashMap::new();

        for (i, gene) in genes.iter().enumerate() {
            let chrom = normalize_chrom(&gene.chrom).into_owned();
            let span = (
                gene.bin_map.offset,
                gene.bin_map.offset + gene.bin_map.bins.len() as u64,
            );
            let entry = chrom_spans.entry(chrom.clone()).or_insert((u64::MAX, 0));
            entry.0 = entry.0.min(span.0);
            entry.1 = entry.1.max(span.1);

            let current_max = chrom_max.entry(chrom.clone()).or_insert(0);
            *current_max = (*current_max).max(gene.representative.genomic_span().1 + 10000);

            genes_by_chrom.entry(chrom).or_default().push(i);
        }

        let mut feature_sizes = HashMap::new();
        if calculate_features {
            // Pre-calculate feature sizes for the whole annotation
            let maps = build_feature_maps(&genes, &chrom_max);
            for map in maps.values() {
                for win in &map.data {
                    for &val in win {
                        let region = match val {
                            0 => RegionType::CdsExon,
                            1 => RegionType::Utr5Exon,
                            2 => RegionType::Utr3Exon,
                            3 => RegionType::Intron,
                            4 => RegionType::TssUp1kb,
                            5 => RegionType::TssUp5kb,
                            6 => RegionType::TssUp10kb,
                            7 => RegionType::TesDown1kb,
                            8 => RegionType::TesDown5kb,
                            9 => RegionType::TesDown10kb,
                            _ => RegionType::Intergenic,
                        };
                        *feature_sizes.entry(region).or_insert(0) += 1;
                    }
                }
            }
        }

        Self {
            version: 3,
            genes,
            feature_sizes,
            chrom_spans,
            genes_by_chrom,
        }
    }

    pub fn build_dense_map(&self, chrom: &str) -> Option<DenseMap> {
        let chrom_norm = normalize_chrom(chrom);
        let (start, end) = self.chrom_spans.get(chrom_norm.as_ref())?;
        self.build_dense_map_for_range(chrom, *start, *end)
    }

    pub fn build_dense_map_for_range(
        &self,
        chrom: &str,
        range_start: u64,
        range_end: u64,
    ) -> Option<DenseMap> {
        let chrom_norm = normalize_chrom(chrom);
        let gene_indices = self.genes_by_chrom.get(chrom_norm.as_ref())?;
        let len = (range_end - range_start) as usize;

        // Stage 1: Single array to track hits.
        let mut map = vec![0u32; len];
        let mut multi_map_temp: HashMap<usize, Vec<u32>> = HashMap::new();

        for &i in gene_indices {
            let gene = &self.genes[i];
            let g_idx = i as u32;
            let g_offset = gene.bin_map.offset;
            let g_len = gene.bin_map.bins.len() as u64;

            // Fast skip if gene is entirely outside the requested range
            if g_offset + g_len <= range_start || g_offset >= range_end {
                continue;
            }

            for (local_pos, &s_pos) in gene.bin_map.bins.iter().enumerate() {
                if s_pos != u32::MAX {
                    let global_pos = g_offset + local_pos as u64;
                    if global_pos >= range_start && global_pos < range_end {
                        let idx = (global_pos - range_start) as usize;
                        let current = map[idx];
                        if current == 0 {
                            map[idx] = g_idx + 1;
                        } else if current == u32::MAX {
                            multi_map_temp.get_mut(&idx).unwrap().push(g_idx);
                        } else {
                            let prev_g_idx = current - 1;
                            map[idx] = u32::MAX;
                            multi_map_temp.insert(idx, vec![prev_g_idx, g_idx]);
                        }
                    }
                }
            }
        }

        // Stage 2: Finalize into the compact DenseMap format
        let mut multi_data = Vec::new();
        let mut multi_offsets = Vec::new();
        multi_offsets.push(0u32);

        // Sort positions to ensure stable multi_map indexing and potentially better locality
        let mut multi_positions: Vec<_> = multi_map_temp.keys().cloned().collect();
        multi_positions.sort_unstable();

        for idx in multi_positions {
            let hits = multi_map_temp.get(&idx).unwrap();
            let m_idx = (multi_offsets.len() - 1) as u32;
            multi_data.extend_from_slice(hits);
            multi_offsets.push(multi_data.len() as u32);
            map[idx] = m_idx | 0x80000000;
        }

        Some(DenseMap {
            offset: range_start,
            map,
            multi_data,
            multi_offsets,
        })
    }

    pub fn build_feature_map(&self, chrom: &str, size: u64) -> Option<FeatureMap> {
        let chrom_norm = normalize_chrom(chrom);
        let gene_indices = self.genes_by_chrom.get(chrom_norm.as_ref())?;
        let mut map = FeatureMap::new(chrom.to_string(), size, 1_000_000);

        for &i in gene_indices {
            crate::analysis::read_distribution::apply_gene_to_map(&self.genes[i], &mut map);
        }

        Some(map)
    }

    pub fn reset_coverage(&self) {
        for gene in self.genes.iter() {
            gene.reset_coverage();
        }
    }

    pub fn save_to_file(&self, path: &str) -> anyhow::Result<()> {
        let f = std::fs::File::create(path)?;
        let writer = std::io::BufWriter::new(f);
        bincode::serialize_into(writer, self)?;
        Ok(())
    }

    pub fn load_from_file(path: &str, max_3p_dist: usize, bin_size: usize) -> anyhow::Result<Self> {
        let f = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(f);
        let mut index: Self = bincode::deserialize_from(reader)?;

        if index.version < 3 {
            anyhow::bail!(
                "Incompatible index version ({} < 3). Please delete old index and regenerate.",
                index.version
            );
        }

        let genes = Arc::get_mut(&mut index.genes).ok_or_else(|| anyhow::anyhow!("Arc busy"))?;
        for gene in genes {
            gene.reinitialize_atomics(max_3p_dist, bin_size);
        }

        Ok(index)
    }
    pub fn estimate_dense_map_memory(&self) -> u64 {
        let mut total_bytes = 0;
        for (start, end) in self.chrom_spans.values() {
            let len = (end - start) as u64;
            // map_base (4 bytes) + final map (4 bytes).
            total_bytes += len * 8;
        }
        total_bytes
    }
}

impl DenseMap {
    #[inline(always)]
    pub fn get_hits(&self, pos: u64) -> Hits<'_> {
        if pos < self.offset {
            return Hits::None;
        }
        let idx = (pos - self.offset) as usize;
        if idx >= self.map.len() {
            return Hits::None;
        }

        let val = self.map[idx];
        if val == 0 {
            Hits::None
        } else if val < 0x80000000 {
            Hits::Single(val - 1)
        } else {
            let m_idx = (val & 0x7FFFFFFF) as usize;
            let start = self.multi_offsets[m_idx] as usize;
            let end = self.multi_offsets[m_idx + 1] as usize;
            Hits::Multi(&self.multi_data[start..end])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::Gene;
    use crate::models::gene::GeneBinMap;

    #[test]
    fn test_build_dense_map_for_range() {
        use crate::models::Transcript;
        let mut gene1 = Gene {
            id: "g1".to_string(),
            name: None,
            chrom: "chr1".to_string(),
            biotype: None,
            total_len: 5,
            representative: Transcript {
                chrom: "chr1".to_string(),
                ..Transcript::default()
            },
            bin_map: Arc::new(GeneBinMap {
                offset: 100,
                bins: vec![0, 1, 2, 3, 4],
            }),
            counts_3p: Arc::new(Vec::new()),
            counts_percentile: Arc::new(Vec::new()),
        };

        let mut gene2 = Gene {
            id: "g2".to_string(),
            name: None,
            chrom: "chr1".to_string(),
            biotype: None,
            total_len: 5,
            representative: Transcript {
                chrom: "chr1".to_string(),
                ..Transcript::default()
            },
            bin_map: Arc::new(GeneBinMap {
                offset: 103,
                bins: vec![0, 1, 2, 3, 4],
            }),
            counts_3p: Arc::new(Vec::new()),
            counts_percentile: Arc::new(Vec::new()),
        };

        let index = AnnotationIndex::new(vec![gene1, gene2], false);

        // Test range fully containing both
        let dense = index.build_dense_map_for_range("chr1", 90, 110).unwrap();
        assert_eq!(dense.offset, 90);
        assert_eq!(dense.map.len(), 20);

        // pos 100 (idx 10): gene1 only
        match dense.get_hits(100) {
            Hits::Single(0) => {}
            _ => panic!("Expected single hit for gene1 at 100"),
        }

        // pos 104 (idx 14): gene1 and gene2
        match dense.get_hits(104) {
            Hits::Multi(hits) => {
                assert_eq!(hits.len(), 2);
                assert!(hits.contains(&0));
                assert!(hits.contains(&1));
            }
            _ => panic!("Expected multi hit at 104"),
        }

        // pos 107 (idx 17): gene2 only
        match dense.get_hits(107) {
            Hits::Single(1) => {}
            _ => panic!("Expected single hit for gene2 at 107"),
        }

        // Test range clipping
        let dense_clipped = index.build_dense_map_for_range("chr1", 104, 106).unwrap();
        assert_eq!(dense_clipped.offset, 104);
        assert_eq!(dense_clipped.map.len(), 2);
        
        // pos 104 in clipped map (idx 0)
        match dense_clipped.get_hits(104) {
            Hits::Multi(_) => {}
            _ => panic!("Expected multi hit at 104 in clipped map"),
        }
    }
}
