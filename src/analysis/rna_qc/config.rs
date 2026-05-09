use crate::analysis::types::AnalysisType;
use crate::io::annotation::IsoformSelect;
use clap::ValueEnum;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, Default, Serialize, Deserialize)]
pub enum DenseMapScope {
    All,
    Chrom,
    #[default]
    Chunk,
}

#[derive(Clone, Debug, Default)]
pub struct RnaQcConfig {
    pub input: Vec<String>,
    pub annotation: String,
    pub annotation_format: String,
    pub output: String,
    pub threads: usize,
    pub analysis: Vec<AnalysisType>,
    pub three_prime_cluster_window: u64,
    pub min_length: u64,
    pub normalization_bp: usize,
    pub gene_id_delimiter: Option<char>,
    pub gene_id_regex: Option<String>,
    pub mapq: u8,
    pub max_3p_dist: usize,
    pub min_support: usize,
    pub no_plot: bool,
    pub snap_qc: bool,
    pub snap_genes: Vec<String>,
    pub snap_flank: u64,
    pub snap_max_reads: usize,
    pub r2_only: bool,
    pub biotype: String,
    pub save_index: bool,
    pub load_index: bool,
    pub ends: bool,
    pub transcript_centric: bool,
    pub plus: bool,
    pub isoform_select: IsoformSelect,
    pub coverage_weighted: bool,
    pub strict_cluster: bool,
    pub stratify_length: bool,
    pub qc_sample_size: usize,
    pub rdna_bed: Option<String>,
    pub rdna_contigs: Vec<String>,
    pub step_size: usize,
    pub three_prime_bin_size: usize,
    pub three_prime_min_anchor_count: u64,
    pub three_prime_min_anchor_mean: f64,
    pub three_prime_min_anchor_nonzero_bins: usize,
    pub three_prime_max_ratio: f64,
    pub write_counts: bool,
    pub dense_map_workers: usize,
    pub dense_map_scope: DenseMapScope,
    pub dense_map_chunk_size: u64,
}
