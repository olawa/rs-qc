pub mod fragment;
pub mod gene;
pub mod transcript;

pub use gene::Gene;
pub use transcript::{Exon, Transcript};

/// Normalizes chromosome names by stripping 'chr' prefix and converting to lowercase.
/// This ensures consistency across different annotation and BAM sources.
pub fn normalize_chrom(name: &str) -> String {
    let n = name.to_lowercase();
    if n.starts_with("chr") {
        n[3..].to_string()
    } else {
        n
    }
}
