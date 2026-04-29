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
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DenseMap {
    pub offset: u64,
    pub map: Vec<u32>, // 0 = None, N < 0x80000000 = genes[N-1] + 1, N >= 0x80000000 = multi_map index
    pub multi_map: Vec<Vec<u32>>,
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

        for gene in genes.iter() {
            let chrom = normalize_chrom(&gene.chrom);
            let span = (
                gene.bin_map.offset,
                gene.bin_map.offset + gene.bin_map.bins.len() as u64,
            );
            let entry = chrom_spans.entry(chrom.clone()).or_insert((u64::MAX, 0));
            entry.0 = entry.0.min(span.0);
            entry.1 = entry.1.max(span.1);

            let current_max = chrom_max.entry(chrom).or_insert(0);
            *current_max = (*current_max).max(gene.representative.genomic_span().1 + 10000);
            // include TES buffer
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
            version: 2,
            genes,
            feature_sizes,
            chrom_spans,
        }
    }

    pub fn build_dense_map(&self, chrom: &str) -> Option<DenseMap> {
        let chrom_norm = normalize_chrom(chrom);
        let (start, end) = self.chrom_spans.get(&chrom_norm)?;
        let len = (end - start) as usize;

        // Stage 1: Single array to track hits.
        // 0: no hit
        // g_idx + 1: single hit
        // u32::MAX: multiple hits (indices stored in multi_map_temp)
        let mut map_base = vec![0u32; len];
        let mut multi_map_temp: HashMap<usize, Vec<u32>> = HashMap::new();

        for (i, gene) in self.genes.iter().enumerate() {
            if normalize_chrom(&gene.chrom) != chrom_norm {
                continue;
            }

            let g_idx = i as u32;
            let g_offset = gene.bin_map.offset;

            for (local_pos, &s_pos) in gene.bin_map.bins.iter().enumerate() {
                if s_pos != u32::MAX {
                    let global_pos = g_offset + local_pos as u64;
                    if global_pos >= *start && global_pos < *end {
                        let idx = (global_pos - *start) as usize;
                        let current = map_base[idx];
                        if current == 0 {
                            map_base[idx] = g_idx + 1;
                        } else if current == u32::MAX {
                            multi_map_temp.get_mut(&idx).unwrap().push(g_idx);
                        } else {
                            // First time an overlap is detected for this base
                            let prev_g_idx = current - 1;
                            map_base[idx] = u32::MAX;
                            multi_map_temp.insert(idx, vec![prev_g_idx, g_idx]);
                        }
                    }
                }
            }
        }

        // Stage 2: Finalize into the compact DenseMap format
        let mut multi_map = Vec::new();
        let mut map = vec![0u32; len];

        // We need a stable mapping from position to multi_map index
        let mut pos_to_m_idx: HashMap<usize, u32> = HashMap::new();

        for (idx, hits) in multi_map_temp {
            let m_idx = multi_map.len() as u32;
            multi_map.push(hits);
            pos_to_m_idx.insert(idx, m_idx | 0x80000000);
        }

        for (idx, val) in map_base.iter().enumerate() {
            if *val == u32::MAX {
                map[idx] = *pos_to_m_idx.get(&idx).unwrap();
            } else {
                map[idx] = *val;
            }
        }

        Some(DenseMap {
            offset: *start,
            map,
            multi_map,
        })
    }

    pub fn build_feature_map(&self, chrom: &str, size: u64) -> Option<FeatureMap> {
        let chrom_norm = normalize_chrom(chrom);
        let mut map = FeatureMap::new(chrom.to_string(), size, 1_000_000);

        for gene in self.genes.iter() {
            if normalize_chrom(&gene.chrom) == chrom_norm {
                crate::analysis::read_distribution::apply_gene_to_map(gene, &mut map);
            }
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

    pub fn load_from_file(path: &str, max_3p_dist: usize) -> anyhow::Result<Self> {
        let f = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(f);
        let mut index: Self = bincode::deserialize_from(reader)?;

        if index.version < 2 {
            anyhow::bail!(
                "Incompatible index version ({} < 2). Please delete old index and regenerate.",
                index.version
            );
        }

        let genes = Arc::get_mut(&mut index.genes).ok_or_else(|| anyhow::anyhow!("Arc busy"))?;
        for gene in genes {
            gene.reinitialize_atomics(max_3p_dist);
        }

        Ok(index)
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
            Hits::Multi(&self.multi_map[m_idx])
        }
    }
}
