pub mod fragment;
pub mod gene;
pub mod transcript;

pub use gene::Gene;
pub use transcript::{Exon, Transcript};

use std::borrow::Cow;

/// Normalizes chromosome names by stripping 'chr' prefix and converting to lowercase.
/// This ensures consistency across different annotation and BAM sources.
pub fn normalize_chrom(name: &str) -> Cow<'_, str> {
    if name.starts_with("chr") || name.chars().any(|c| c.is_uppercase()) {
        let n = name.to_lowercase();
        if n.starts_with("chr") {
            Cow::Owned(n[3..].to_string())
        } else {
            Cow::Owned(n)
        }
    } else {
        Cow::Borrowed(name)
    }
}

pub fn is_autosomal_chrom(name: &str) -> bool {
    normalize_chrom(name)
        .parse::<u8>()
        .map(|n| (1..=22).contains(&n))
        .unwrap_or(false)
}
