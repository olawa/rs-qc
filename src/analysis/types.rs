use clap::ValueEnum;
use std::fmt;

#[derive(ValueEnum, Clone, Debug, PartialEq, Eq, Copy, Hash)]
pub enum AnalysisType {
    #[value(alias = "genebody_coverage", alias = "genebody")]
    GeneBody,
    #[value(alias = "3prime", alias = "three_prime")]
    ThreePrime,
    #[value(alias = "read_distribution", alias = "distribution")]
    Distribution,
    #[value(alias = "rna_qc", alias = "qc")]
    Qc,
}

impl fmt::Display for AnalysisType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GeneBody => write!(f, "genebody"),
            Self::ThreePrime => write!(f, "three_prime"),
            Self::Distribution => write!(f, "distribution"),
            Self::Qc => write!(f, "qc"),
        }
    }
}
