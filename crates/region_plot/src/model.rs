#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct RegionPlot {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    /// Reference/consensus bases for the viewed interval. This is used for the
    /// rsnap-style colored top reference strip.
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
    /// Rich per-base pileup. Prefer this for rsnap-style plots.
    #[cfg_attr(feature = "serde", serde(default))]
    pub pileup: Vec<BasePileup>,
    /// Legacy/minimal coverage input. Used when `pileup` is empty.
    #[cfg_attr(feature = "serde", serde(default))]
    pub coverage: Vec<CoveragePoint>,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct CoveragePoint {
    pub pos: i64,
    pub depth: u32,
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BasePileup {
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
    pub del: u32,
    pub ins: u32,
    pub ins_len: u32,
    pub total: u32,
    pub del_spanning: u32,
    pub del_starts: u32,
}

impl BasePileup {
    pub fn depth(&self) -> u32 {
        self.total.max(self.a + self.c + self.g + self.t + self.n)
    }

    pub fn base_count(&self, base: u8) -> u32 {
        match base.to_ascii_uppercase() {
            b'A' => self.a,
            b'C' => self.c,
            b'G' => self.g,
            b'T' => self.t,
            b'N' => self.n,
            _ => 0,
        }
    }
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
    Mismatch {
        ref_start: i64,
        len: i64,
        query_start: usize,
        base: u8,
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
    Skip {
        ref_start: i64,
        len: i64,
    },
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct BaseModification {
    pub ref_pos: i64,
    pub code: String,
    pub probability: Option<f32>,
}
