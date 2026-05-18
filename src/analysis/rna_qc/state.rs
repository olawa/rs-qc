use crate::analysis::feature_index::FeatureIndex;
use crate::analysis::qc::InlineQcState;
use anyhow::Result;
use std::fs::File;
use std::io::Write;

#[derive(Debug, Clone)]
pub(crate) struct RnaWorkerState {
    pub(crate) records_seen: u64,
    pub(crate) fail_unmapped: u64,
    pub(crate) fail_secondary: u64,
    pub(crate) fail_qc: u64,
    pub(crate) fail_mapq: u64,
    pub(crate) overlaps_found: u64,
    pub(crate) total_tags: u64,
    pub(crate) aligned_qc_reads: u64,
    pub(crate) mtdna_reads: u64,
    pub(crate) rdna_reads: u64,
    pub(crate) read_dist_counts: [u64; 12],
    pub(crate) unknown_chrom_reads: u64,
    pub(crate) qc: InlineQcState,
    pub(crate) config: super::config::RnaQcConfig,
    pub(crate) distribution_transcripts_loaded: usize,
    pub(crate) distribution_genes_loaded: usize,
    pub(crate) distribution_biotypes_included: String,
    pub(crate) distribution_total_classified_reads: u64,
    pub(crate) distribution_exonic_reads: u64,
    pub(crate) distribution_intronic_reads: u64,
    pub(crate) distribution_flank_reads: u64,
    pub(crate) distribution_intergenic_reads: u64,
    pub(crate) splice_junctions: std::collections::HashMap<crate::analysis::splice_junction::JunctionKey, crate::analysis::splice_junction::ObservedJunction>,
}

impl RnaWorkerState {
    pub(crate) fn new(config: &super::config::RnaQcConfig, qc_sample_size: usize) -> Self {
        Self {
            records_seen: 0,
            fail_unmapped: 0,
            fail_secondary: 0,
            fail_qc: 0,
            fail_mapq: 0,
            total_tags: 0,
            overlaps_found: 0,
            aligned_qc_reads: 0,
            mtdna_reads: 0,
            rdna_reads: 0,
            read_dist_counts: [0; 12],
            unknown_chrom_reads: 0,
            qc: InlineQcState::new(qc_sample_size),
            config: config.clone(),
            distribution_transcripts_loaded: 0,
            distribution_genes_loaded: 0,
            distribution_biotypes_included: String::new(),
            distribution_total_classified_reads: 0,
            distribution_exonic_reads: 0,
            distribution_intronic_reads: 0,
            distribution_flank_reads: 0,
            distribution_intergenic_reads: 0,
            splice_junctions: std::collections::HashMap::new(),
        }
    }

    pub(crate) fn merge(mut self, other: Self) -> Self {
        self.records_seen += other.records_seen;
        self.fail_unmapped += other.fail_unmapped;
        self.fail_secondary += other.fail_secondary;
        self.fail_qc += other.fail_qc;
        self.fail_mapq += other.fail_mapq;
        self.overlaps_found += other.overlaps_found;
        self.total_tags += other.total_tags;
        self.aligned_qc_reads += other.aligned_qc_reads;
        self.mtdna_reads += other.mtdna_reads;
        self.rdna_reads += other.rdna_reads;

        for i in 0..12 {
            self.read_dist_counts[i] += other.read_dist_counts[i];
        }
        self.unknown_chrom_reads += other.unknown_chrom_reads;

        self.distribution_total_classified_reads += other.distribution_total_classified_reads;
        self.distribution_exonic_reads += other.distribution_exonic_reads;
        self.distribution_intronic_reads += other.distribution_intronic_reads;
        self.distribution_flank_reads += other.distribution_flank_reads;
        self.distribution_intergenic_reads += other.distribution_intergenic_reads;

        self.qc.merge_from(other.qc);

        for (key, other_j) in other.splice_junctions {
            self.splice_junctions.entry(key)
                .and_modify(|j| {
                    j.count += other_j.count;
                    j.min_hash = j.min_hash.min(other_j.min_hash);
                })
                .or_insert(other_j);
        }

        self
    }

    pub(crate) fn write_read_distribution_report(
        &self,
        path: &str,
        feature_index: &FeatureIndex,
    ) -> Result<()> {
        let mut f_dist = File::create(path)?;
        use crate::analysis::read_distribution::RegionType;

        writeln!(f_dist, "{:<38}{}", "distribution_transcripts_loaded", self.distribution_transcripts_loaded)?;
        writeln!(f_dist, "{:<38}{}", "distribution_genes_loaded", self.distribution_genes_loaded)?;
        writeln!(f_dist, "{:<38}{}", "distribution_biotypes_included", self.distribution_biotypes_included)?;
        writeln!(f_dist, "{:<38}{}", "distribution_total_classified_reads", self.distribution_total_classified_reads)?;
        writeln!(f_dist, "{:<38}{}", "distribution_exonic_reads", self.distribution_exonic_reads)?;
        writeln!(f_dist, "{:<38}{}", "distribution_intronic_reads", self.distribution_intronic_reads)?;
        writeln!(f_dist, "{:<38}{}", "distribution_flank_reads", self.distribution_flank_reads)?;
        writeln!(f_dist, "{:<38}{}", "distribution_intergenic_reads", self.distribution_intergenic_reads)?;
        writeln!(
            f_dist,
            "====================================================================="
        )?;

        let total_assigned: u64 = (0..12)
            .filter(|&i| i != RegionType::Intergenic as usize)
            .map(|i| self.read_dist_counts[i])
            .sum();

        let passed_filters = self.records_seen
            - self.fail_unmapped
            - self.fail_secondary
            - self.fail_qc
            - self.fail_mapq;

        writeln!(f_dist, "{:<30}{}", "Total Reads", passed_filters)?;
        writeln!(f_dist, "{:<30}{}", "Total Tags", self.total_tags)?;
        writeln!(f_dist, "{:<30}{}", "Total Assigned Tags", total_assigned)?;
        writeln!(
            f_dist,
            "====================================================================="
        )?;
        writeln!(
            f_dist,
            "{:<20}{:<20}{:<20}{:<20}",
            "Group", "Total_bases", "Tag_count", "Tags/Kb"
        )?;

        let groups = [
            RegionType::CdsExon,
            RegionType::Utr5Exon,
            RegionType::Utr3Exon,
            RegionType::Exon,
            RegionType::Intron,
            RegionType::TssUp1kb,
            RegionType::TssUp5kb,
            RegionType::TssUp10kb,
            RegionType::TesDown1kb,
            RegionType::TesDown5kb,
            RegionType::TesDown10kb,
            RegionType::Intergenic,
        ];

        for group in groups {
            let size = feature_index
                .feature_sizes
                .get(&group)
                .cloned()
                .unwrap_or(0);
            let count = self.read_dist_counts[group as usize];
            let tags_per_kb = if size > 0 {
                (count as f64 * 1000.0) / (size as f64)
            } else {
                0.0
            };
            writeln!(
                f_dist,
                "{:<20}{:<20}{:<20}{:<18.2}",
                group.to_string(),
                size,
                count,
                tags_per_kb
            )?;
        }
        writeln!(
            f_dist,
            "====================================================================="
        )?;
        Ok(())
    }
}
