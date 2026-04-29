# rs-qc Quick Start

This page gives a short command-by-command tour of the current CLI.

## FASTQ

Streaming FASTQ QC with per-base quality, GC, adapters, duplication estimates, and k-mers.

```bash
rs-qc fastq -i reads_R1.fastq.gz reads_R2.fastq.gz --paired -o sample_fastq
```

Useful options:

- `--sample-size` for sampled duplication and k-mer tracking
- `--kmer-size` to tune the k-mer screen
- `--top-n` to control how many overrepresented sequences are reported
- `--no-kmers` to skip k-mer counting

## Alignment

General BAM/CRAM QC with flag counts, MAPQ, insert size, clipping, contig summaries, and `de:f` accuracy.

```bash
rs-qc align -i sample.bam -o sample_align -t 8
```

Useful options:

- `--mapq` to set the filtering threshold used in summaries
- `--threads` to control BGZF worker count

## RNA

RNA-seq QC with gene body coverage, strandness, inner distance, read distribution, and contamination-aware summaries.

```bash
rs-qc rna \
  -i sample.bam \
  -a annotation.gtf \
  -o sample_rna \
  --analysis gene-body,three-prime,distribution,qc
```

Useful options:

- `--ends` for 3' end-focused analysis
- `--qc-sample-size` for strandedness and inner-distance sampling
- `--rdna-bed` for explicit rDNA interval annotation
- `--rdna-contigs` for contig-name-based rDNA detection

## Report

Render a unified JSON and HTML report from existing module summaries.

```bash
rs-qc report -i sample_fastq sample_align sample_rna -o sample_report
```

Useful inputs:

- an output prefix like `sample_fastq`
- a summary JSON path like `sample.fastq.summary.json`
- a mix of module prefixes and summary JSON files

## Planned Modules

These are exposed in the CLI shell already, but are not implemented yet:

- `rs-qc dna`
- `rs-qc atac`
- `rs-qc contam`

## Output

All modules write plain-text or TSV outputs using the output prefix you choose.

That design is deliberate: the long-term plan is to render JSON and HTML reports from the same structured metrics without re-scanning the inputs.
