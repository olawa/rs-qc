use crate::models::Gene;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
pub enum RegionType {
    CdsExon = 0,
    Utr5Exon = 1,
    Utr3Exon = 2,
    Intron = 3,
    TssUp1kb = 4,
    TssUp5kb = 5,
    TssUp10kb = 6,
    TesDown1kb = 7,
    TesDown5kb = 8,
    TesDown10kb = 9,
    Intergenic = 10,
    Exon = 11,
}

impl std::fmt::Display for RegionType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CdsExon => write!(f, "CDS_Exons"),
            Self::Utr5Exon => write!(f, "5'UTR_Exons"),
            Self::Utr3Exon => write!(f, "3'UTR_Exons"),
            Self::Intron => write!(f, "Introns"),
            Self::TssUp1kb => write!(f, "TSS_up_1kb"),
            Self::TssUp5kb => write!(f, "TSS_up_5kb"),
            Self::TssUp10kb => write!(f, "TSS_up_10kb"),
            Self::TesDown1kb => write!(f, "TES_down_1kb"),
            Self::TesDown5kb => write!(f, "TES_down_5kb"),
            Self::TesDown10kb => write!(f, "TES_down_10kb"),
            Self::Intergenic => write!(f, "Intergenic"),
            Self::Exon => write!(f, "Exons"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureMap {
    pub chrom: String,
    pub window_size: u64,
    pub data: Vec<Vec<u8>>,
}

impl FeatureMap {
    pub fn new(chrom: String, size: u64, window_size: u64) -> Self {
        let n_windows = (size as f64 / window_size as f64).ceil() as usize;
        Self {
            chrom,
            window_size,
            data: vec![vec![RegionType::Intergenic as u8; window_size as usize]; n_windows],
        }
    }

    pub fn set_region(&mut self, start: u64, end: u64, region: RegionType) {
        let region_val = region as u8;
        for pos in start..end {
            let win_idx = (pos / self.window_size) as usize;
            let offset = (pos % self.window_size) as usize;
            if win_idx < self.data.len() {
                let current = self.data[win_idx][offset];
                if region_val < current {
                    self.data[win_idx][offset] = region_val;
                }
            }
        }
    }

    pub fn get_region(&self, pos: u64) -> RegionType {
        let win_idx = (pos / self.window_size) as usize;
        let offset = (pos % self.window_size) as usize;
        if win_idx < self.data.len() && offset < self.data[win_idx].len() {
            match self.data[win_idx][offset] {
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
            }
        } else {
            RegionType::Intergenic
        }
    }
}

pub fn apply_gene_to_map(gene: &Gene, map: &mut FeatureMap) {
    if crate::models::normalize_chrom(&gene.chrom) != map.chrom {
        return;
    }
    // ... rest of function remains same
    let tx = &gene.representative;

    for exon in &tx.exons {
        if let (Some(cds_s), Some(cds_e)) = (tx.cds_start, tx.cds_end) {
            let overlap_s = exon.start.max(cds_s);
            let overlap_e = exon.end.min(cds_e);
            if overlap_s < overlap_e {
                map.set_region(overlap_s, overlap_e, RegionType::CdsExon);
            }

            if tx.strand == '+' {
                let utr5_e = exon.end.min(cds_s);
                if exon.start < utr5_e {
                    map.set_region(exon.start, utr5_e, RegionType::Utr5Exon);
                }

                let utr3_s = exon.start.max(cds_e);
                if utr3_s < exon.end {
                    map.set_region(utr3_s, exon.end, RegionType::Utr3Exon);
                }
            } else {
                let utr5_s = exon.start.max(cds_e);
                if utr5_s < exon.end {
                    map.set_region(utr5_s, exon.end, RegionType::Utr5Exon);
                }

                let utr3_e = exon.end.min(cds_s);
                if exon.start < utr3_e {
                    map.set_region(exon.start, utr3_e, RegionType::Utr3Exon);
                }
            }
        } else {
            map.set_region(exon.start, exon.end, RegionType::Utr5Exon);
        }
    }

    for i in 0..tx.exons.len().saturating_sub(1) {
        let intron_s = tx.exons[i].end;
        let intron_e = tx.exons[i + 1].start;
        if intron_s < intron_e {
            map.set_region(intron_s, intron_e, RegionType::Intron);
        }
    }

    let tss = tx.five_prime_end();
    let tes = tx.three_prime_end();

    if tx.strand == '+' {
        map.set_region(
            tss.saturating_sub(10000),
            tss.saturating_sub(5000),
            RegionType::TssUp10kb,
        );
        map.set_region(
            tss.saturating_sub(5000),
            tss.saturating_sub(1000),
            RegionType::TssUp5kb,
        );
        map.set_region(tss.saturating_sub(1000), tss, RegionType::TssUp1kb);

        map.set_region(tes, tes + 1000, RegionType::TesDown1kb);
        map.set_region(tes + 1000, tes + 5000, RegionType::TesDown5kb);
        map.set_region(tes + 5000, tes + 10000, RegionType::TesDown10kb);
    } else {
        map.set_region(tss + 5000, tss + 10000, RegionType::TssUp10kb);
        map.set_region(tss + 1000, tss + 5000, RegionType::TssUp5kb);
        map.set_region(tss, tss + 1000, RegionType::TssUp1kb);

        map.set_region(tes.saturating_sub(1000), tes, RegionType::TesDown1kb);
        map.set_region(
            tes.saturating_sub(5000),
            tes.saturating_sub(1000),
            RegionType::TesDown5kb,
        );
        map.set_region(
            tes.saturating_sub(10000),
            tes.saturating_sub(5000),
            RegionType::TesDown10kb,
        );
    }
}

pub fn build_feature_maps(
    genes: &[Gene],
    chrom_sizes: &HashMap<String, u64>,
) -> HashMap<String, FeatureMap> {
    let mut maps = HashMap::new();
    let window_size = 1_000_000;

    for (chrom, &size) in chrom_sizes {
        // chrom is already normalized here because chrom_sizes was built with normalized keys
        maps.insert(
            chrom.clone(),
            FeatureMap::new(chrom.clone(), size, window_size),
        );
    }

    for gene in genes {
        let norm_chrom = crate::models::normalize_chrom(&gene.chrom);
        if let Some(map) = maps.get_mut(&norm_chrom) {
            apply_gene_to_map(gene, map);
        }
    }

    maps
}
