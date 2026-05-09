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
    pub(crate) qc: InlineQcState,
}

impl RnaWorkerState {
    pub(crate) fn new(qc_sample_size: usize) -> Self {
        Self {
            records_seen: 0,
            fail_unmapped: 0,
            fail_secondary: 0,
            fail_qc: 0,
            fail_mapq: 0,
            overlaps_found: 0,
            total_tags: 0,
            aligned_qc_reads: 0,
            mtdna_reads: 0,
            rdna_reads: 0,
            read_dist_counts: [0; 12],
            qc: InlineQcState::new(qc_sample_size),
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

        self.qc.merge_from(other.qc);
        self
    }

    pub(crate) fn write_read_distribution_report(
        &self,
        path: &str,
        feature_index: &FeatureIndex,
    ) -> Result<()> {
        let mut f_dist = File::create(path)?;
        use crate::analysis::read_distribution::RegionType;

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
