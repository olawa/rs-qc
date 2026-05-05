# region_plot

`region_plot` is a small, pure plotting crate for static genomic region reports.
It does not read BAM/CRAM, FASTA, GTF, BED, or VCF files. Upstream tools provide
already-extracted reads, coverage, reference bases, and annotations as plain Rust
data structures.

The crate is intentionally independent from `rsnap` viewer internals. It can be
moved to a separate repository without changing rsnap.

## Example

```rust
use region_plot::{
    render_to_path, BasePileup, CoveragePoint, GeneModel, PlotOptions, RegionPlot, ReadModel,
    ReadSegment, SamplePlotData,
};

let plot = RegionPlot {
    chrom: "chr1".to_string(),
    start: 100,
    end: 220,
    reference: Some(b"ACGTACGTACGTACGTACGTACGTACGT".to_vec()),
    genes: vec![GeneModel {
        name: "GENE1".to_string(),
        start: 120,
        end: 200,
        strand: Some('+'),
        exons: vec![(130, 150), (170, 190)],
    }],
    samples: vec![SamplePlotData {
        name: "sample".to_string(),
        pileup: vec![BasePileup {
            a: 12,
            total: 12,
            ..BasePileup::default()
        }],
        coverage: vec![CoveragePoint { pos: 100, depth: 12 }],
        reads: vec![ReadModel {
            name: "read-1".to_string(),
            start: 110,
            end: 180,
            is_reverse: false,
            mapq: 60,
            segments: vec![ReadSegment::Match {
                ref_start: 110,
                len: 70,
                query_start: 0,
            }],
            bases: None,
            qualities: None,
            haplotype: Some(1),
            modifications: Vec::new(),
        }],
    }],
};

render_to_path(&plot, &PlotOptions::default(), "region.png")?;
# anyhow::Ok(())
```

The preferred coverage input is `SamplePlotData::pileup`, which mirrors the
information rsnap needs for coverage and base tracks. `coverage` is retained as a
minimal fallback for callers that only have depth values.
