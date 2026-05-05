use crate::analysis::index::AnnotationIndex;
use crate::models::{is_autosomal_chrom, normalize_chrom};
use std::collections::HashSet;

pub(crate) fn normalized_name_set(names: &[String]) -> HashSet<String> {
    names
        .iter()
        .map(|n| normalize_chrom(n).into_owned())
        .collect()
}

pub(crate) fn autosomal_coverage_index(index: &AnnotationIndex) -> AnnotationIndex {
    let genes = index
        .genes
        .iter()
        .filter(|gene| is_autosomal_chrom(&gene.chrom))
        .cloned()
        .collect();
    AnnotationIndex::new(genes, false)
}

pub(crate) fn normalize_thread_count(threads: usize) -> usize {
    threads.max(1)
}

pub(crate) fn window_owns_record_start(pos: u64, win_start: u64, win_end: u64) -> bool {
    pos >= win_start && pos < win_end
}

pub(crate) fn is_mtdna_chrom_norm(chrom: &str) -> bool {
    chrom == "m" || chrom == "mt" || chrom == "chrm" || chrom == "chrmt"
}

pub(crate) fn is_rdna_chrom_norm(chrom: &str, rdna_contigs: &HashSet<String>) -> bool {
    rdna_contigs.contains(chrom)
}
