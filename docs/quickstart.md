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
- `--length-bin-size` to tune long-read length histogram bins
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

## DNA

Mosdepth-like coverage summaries over the genome, windows, and optional targets.

```bash
rs-qc dna -i sample.bam -o sample_dna --window-size 100000 --targets panel.bed
```

Useful options:

- `--window-size` for genome windowing
- `--targets` for BED-based target coverage

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
- `--snap-qc` to render marker-gene snapshots for RNA QC
- `--snap-genes` to add custom snapshot genes such as `GAPDH,ACTB,MALAT1`
- `--snap-flank` and `--snap-max-reads` to tune snapshot windows and density

RNA output now includes a compact terminal summary, per-sample `*.rna.qc_summary.svg`,
and optional `*.rna_snapshots.tsv` plus a snapshot directory.

## Snapshot

Static genomic region snapshots with coverage, reads, CIGAR events, optional annotation, and optional reference bases.

```bash
rs-qc snap \
  -i sample.bam \
  -r chr22:28600000-28602000 \
  -a gencode.gtf.gz \
  --reference hg38.fa \
  -o sample.chr22.snap.png
```

Useful options:

- `--max-reads` to cap displayed reads
- `--sample-reads` to reservoir-sample displayed reads when the region is crowded
- `--no-reference` or `--no-genes` to suppress optional tracks
- `--format` to require `png` or `svg` output

## Contamination

Aligned mtDNA/rDNA summaries plus optional exact k-mer screening against rRNA references.

```bash
rs-qc contam \
  -i sample.bam \
  -o sample_contam \
  --rdna-bed rdna.bed \
  --rrna-fasta rrna_refs.fa
```

Useful options:

- `--rdna-contigs` for contig-name-based rDNA/rRNA detection
- `--rrna-fasta` to build an exact rRNA k-mer screen
- `--sample-size` to cap the number of reads used for the k-mer estimate
- `--kmer-scan-all-reads` to screen all reads instead of unmapped/low-MAPQ reads

## Report

Render a unified JSON and HTML report from existing module summaries.

```bash
rs-qc report -i sample_fastq sample_align sample_rna sample_dna -o sample_report
```

Useful inputs:

- an output prefix like `sample_fastq`
- a summary JSON path like `sample.fastq.summary.json`
- a mix of module prefixes and summary JSON files

RNA report sections embed the combined QC summary figure and link any snapshots
that were generated for the sample.

## Planned Modules

These are exposed in the CLI shell already, but are not implemented yet:

- `rs-qc atac`

## Output

All modules write plain-text or TSV outputs using the output prefix you choose.

That design is deliberate: the long-term plan is to render JSON and HTML reports from the same structured metrics without re-scanning the inputs.
