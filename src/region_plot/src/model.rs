#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct RegionPlot {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub reference: Option<Vec<u8>>,
    pub genes: Vec<GeneModel>,
    pub samples: Vec<SamplePlotData>,
}

impl RegionPlot {
    pub fn span(&self) -> i64 {
        (self.end - self.start).max(1)
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct SamplePlotData {
    pub name: String,
    pub reads: Vec<ReadModel>,
    pub coverage: Vec<CoveragePoint>,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct CoveragePoint {
    pub pos: i64,
    pub depth: u32,
}

pub type Strand = char;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct GeneModel {
    pub name: String,
    pub start: i64,
    pub end: i64,
    pub strand: Option<Strand>,
    pub exons: Vec<(i64, i64)>,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct ReadModel {
    pub name: String,
    pub start: i64,
    pub end: i64,
    pub is_reverse: bool,
    pub mapq: u8,
    pub segments: Vec<ReadSegment>,
    pub bases: Option<Vec<u8>>,
    pub qualities: Option<Vec<u8>>,
    pub haplotype: Option<u8>,
    pub modifications: Vec<BaseModification>,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub enum ReadSegment {
    Match {
        ref_start: i64,
        len: i64,
        query_start: usize,
    },
    Ins {
        ref_pos: i64,
        bases: Vec<u8>,
    },
    Del {
        ref_start: i64,
        len: i64,
    },
    SoftClip {
        ref_pos: i64,
        bases: Vec<u8>,
    },
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct BaseModification {
    pub ref_pos: i64,
    pub code: String,
    pub probability: Option<f32>,
}
