# rs-qc

> Rapid sequencing qc for FASTQ, BAM/CRAM, RNA-seq, and more.

- Version: `0.1.0`
- Language: Rust
- Scope: Raw reads, alignment QC, RNA QC, and future assay-specific modules

`rs-qc` stands for **rapid sequencing qc**.

It is a fast Rust toolkit for sequencing quality control across the whole NGS pipeline: raw FASTQ, aligned BAM/CRAM, and assay-specific QC for RNA-seq, DNA-seq, and contamination checks.

The goal is not to be a Rust port of one legacy tool. The goal is to be a single, modern QC binary that stays fast on large files, emits structured metrics, and can grow into a complete report generator.

## Why rs-qc?

- `rapid` because the tools should be fast enough to run during everyday analysis, not only in overnight QC batches.
- `sequencing` because the scope is broader than RNA-seq and should cover the common NGS workflows people actually use.
- `qc` because the output should be practical, assay-aware quality control rather than one-off summary stats.

The project intentionally starts with a strong RNA-seq core, but the shape is broader:

- raw read QC
- general alignment QC
- DNA coverage QC
- RNA QC
- contamination screening
- future report generation

## What It Does Today

`rs-qc` currently includes:

- `rs-qc fastq` for streaming FASTQ QC
- `rs-qc align` for general BAM/CRAM alignment QC
- `rs-qc rna` for RNA-seq QC with RSeQC-style metrics
- `rs-qc dna` for mosdepth-like coverage and breadth QC
- `rs-qc snap` for static genomic region snapshots
- `rs-qc report` for unified JSON and HTML reporting

The CLI also exposes planned subcommands for future expansion:

- `rs-qc atac`
- `rs-qc contam`

## Features

### FASTQ QC

- Per-base quality
- Per-read mean quality
- GC distribution
- Read length distribution
- Per-base A/C/G/T/N composition
- Adapter detection
- Overrepresented sequences
- Duplication estimate from sampled reads
- Poly-A / poly-G tail detection
- K-mer enrichment
- Optional paired-end name sync check

### Alignment QC

- Total, mapped, unmapped, primary, secondary, and supplementary reads
- Duplicate and QC-fail counts
- Proper pair, singleton, discordant, and orphan counts
- MAPQ histogram
- Read-length histogram
- Insert-size histogram
- CIGAR operation summary
- Soft-clipping and hard-clipping summary
- Per-contig read counts
- `de:f` alignment accuracy summaries
- Accuracy histogram with mode, median, mean, and quantiles

### RNA-seq QC

- Gene body coverage
- 3' coverage bias
- Read distribution across CDS, UTR, exon, intron, and intergenic regions
- Strandedness / infer experiment
- Inner distance distribution
- mtDNA fraction
- rDNA fraction
- rRNA / contaminant interval support

## Build

```bash
cargo build --release
```

## Quick Start

### FASTQ

```bash
rs-qc fastq -i reads_R1.fastq.gz reads_R2.fastq.gz --paired -o sample_fastq
```

### Alignment

```bash
rs-qc align -i sample.bam -o sample_align -t 8
```

### RNA-seq

```bash
rs-qc rna \
  -i sample.bam \
  -a annotation.gtf \
  -o sample_rna \
  --analysis gene-body,three-prime,distribution,qc
```

### DNA coverage

```bash
rs-qc dna -i sample.bam -o sample_dna --window-size 100000 --targets panel.bed
```

### Contamination

```bash
rs-qc contam \
  -i sample.bam \
  -o sample_contam \
  --rdna-bed rdna.bed \
  --rrna-fasta rrna_refs.fa \
  --kmer-size 31 \
  --min-kmer-hits 2
```

### Region snapshot

```bash
rs-qc snap \
  -i sample.bam \
  -r chr22:28600000-28602000 \
  -a gencode.gtf.gz \
  --reference hg38.fa \
  -o sample.chr22.snap.png
```

## Output Files

Each module writes plain-text or TSV outputs with the chosen output prefix.

Examples:

- `sample.fastq.summary.txt`
- `sample.fastq.summary.json`
- `sample.fastq.per_base.tsv`
- `sample.fastq.length_distribution.tsv`
- `sample.fastq.length_bins.tsv`
- `sample.fastq.length_distribution.svg`
- `sample.align.summary.txt`
- `sample.align.summary.json`
- `sample.align.mapq.tsv`
- `sample.rna_qc.txt`
- `sample.rna.summary.json`
- `sample.inner_distance.tsv`
- `sample.geneBodyCoverage.txt`
- `sample.summary.json`
- `sample.report.html`
- `sample.dna.summary.tsv`
- `sample.dna.summary.json`
- `sample.dna.depth_hist.tsv`
- `sample.dna.contigs.tsv`
- `sample.dna.windows.tsv`
- `sample.dna.targets.tsv`
- `sample.contam.summary.tsv`
- `sample.contam.summary.json`
- `sample.chr22.snap.png`

This structured output is intentional: it makes the future JSON and HTML report layer much easier to build without re-running the scan.

## Design Notes

- `rs-qc` is built for large files and tries to keep the hot path simple.
- Metrics are collected into typed structs first, then rendered to text/TSV.
- The code is organized so shared BAM scanning can be reused across modules instead of re-implementing progress bars and indexing logic everywhere.
- The project currently favors deterministic, streaming QC over fancy post-hoc inference.

## Roadmap

Planned next steps include:

- `rs-qc atac` for ATAC/ChIP-style metrics
- species-level contamination sketches or marker-kmer screens
- richer HTML report cards and warning badges

## Notes

- Input FASTQ files can be gzipped.
- BAM input is currently the most mature alignment path.
- CRAM support is planned, but BAM is the current focus.
- DNA coverage uses a simple coordinate-sweep accumulator, not the RNA transcript index.
- Region snapshots keep `region_plot` as a pure renderer; BAM, FASTA, and annotation extraction live in `rs-qc`.

## More Docs

- [Quick Start](docs/quickstart.md)
