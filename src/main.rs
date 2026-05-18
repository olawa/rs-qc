mod analysis;
mod io;
mod models;
mod stats;

use crate::analysis::alignment_qc::{
    run_alignment_qc, sample_name_from_alignment_path, AlignmentQcConfig,
};
use crate::analysis::contamination::{run_contamination_qc, ContaminationQcConfig};
use crate::analysis::dna_qc::{run_dna_qc, DnaQcConfig};
use crate::analysis::fastq_qc::{run_fastq_qc, sample_name_from_path, FastqQcConfig};
use crate::analysis::report::{build_document, resolve_input_files, write_report};
use crate::analysis::rna_qc::{run_rna, RnaQcConfig};
use crate::analysis::snapshot::{
    resolve_snapshot_region, run_snapshot, validate_output_format, SnapshotConfig,
    SnapshotOutputFormat,
};
use crate::analysis::types::AnalysisType;
use crate::io::annotation::AnnotationFormat;

use anyhow::{bail, Result};
use clap::{Args as ClapArgs, Parser, Subcommand, ValueEnum};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Rapid sequencing QC for FASTQ, BAM/CRAM, and assay-specific NGS metrics"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// RNA-seq QC, including gene-body coverage and RSeQC-like annotation metrics.
    Rna(RnaArgs),
    /// FASTQ QC metrics such as quality, GC, adapters, duplication, and k-mers.
    Fastq(FastqArgs),
    /// General alignment QC for BAM/CRAM files.
    Align(AlignArgs),
    /// DNA coverage QC with mosdepth-like depth and breadth metrics.
    Dna(DnaArgs),
    /// ATAC/ChIP/cfDNA-style assay metrics.
    Atac(PlannedArgs),
    /// Contamination metrics from aligned intervals, contigs, and optional k-mer screens.
    Contam(ContamArgs),
    /// Render JSON/TSV/SVG metrics into a unified report.
    Report(ReportArgs),
    /// Render a static genomic region snapshot from BAM, annotation, and reference data.
    Snap(SnapArgs),
}

#[derive(ClapArgs, Debug)]
struct PlannedArgs {
    /// Output prefix for future module outputs.
    #[arg(short, long, default_value = "rs-qc")]
    output: String,
}

#[derive(ClapArgs, Debug)]
struct ContamArgs {
    #[arg(short, long, num_args = 1.., required = true)]
    input: Vec<String>,
    #[arg(short, long)]
    output: Option<String>,
    #[arg(long)]
    rdna_bed: Option<String>,
    #[arg(long, value_delimiter = ',', default_values_t = crate::analysis::rna_qc::default_rdna_contigs())]
    rdna_contigs: Vec<String>,
    #[arg(long)]
    rrna_fasta: Option<String>,
    #[arg(long, default_value_t = 31)]
    kmer_size: usize,
    #[arg(long, default_value_t = 2)]
    min_kmer_hits: usize,
    #[arg(long, default_value_t = 1_000_000)]
    sample_size: usize,
    #[arg(long, default_value_t = 10)]
    mapq: u8,
    #[arg(long, default_value_t = false)]
    kmer_scan_all_reads: bool,
}

#[derive(ClapArgs, Debug)]
struct DnaArgs {
    #[arg(short, long, num_args = 1..)]
    input: Vec<String>,
    #[arg(short, long)]
    output: Option<String>,
    #[arg(long, default_value_t = 10)]
    mapq: u8,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
    #[arg(long, default_value_t = 100_000)]
    window_size: u32,
    #[arg(long)]
    targets: Option<String>,
    #[arg(long, value_delimiter = ',', default_values_t = vec![1, 5, 10, 20, 30])]
    thresholds: Vec<u32>,
    #[arg(long, default_value_t = 10)]
    callable_depth: u32,
    #[arg(long, default_value_t = false)]
    include_duplicates: bool,
    #[arg(short, long)]
    reference: Option<String>,
    #[arg(long)]
    annotation: Option<String>,
    #[arg(long, default_value_t = 5)]
    low_cov_threshold: u32,
    #[arg(long, default_value_t = 0)]
    snap_lowcov: u32,
}

#[derive(ClapArgs, Debug)]
struct ReportArgs {
    #[arg(short, long, num_args = 1.., required = true)]
    input: Vec<String>,
    #[arg(short, long)]
    output: Option<String>,
}

#[derive(ValueEnum, Clone, Debug, PartialEq, Eq, Default)]
pub enum SnapAnnotationFormat {
    #[default]
    Auto,
    Gtf,
    Bed12,
}

#[derive(ValueEnum, Clone, Debug, PartialEq, Eq, Default)]
pub enum SnapFormat {
    #[default]
    Auto,
    Png,
    Svg,
}

#[derive(ValueEnum, Clone, Debug, PartialEq, Eq, Default)]
pub enum SnapBaseTrackMode {
    #[default]
    Top,
    Bottom,
    Both,
    None,
}

#[derive(ClapArgs, Debug)]
struct SnapArgs {
    #[arg(short = 'i', long = "bam")]
    bam: String,
    #[arg(short, long)]
    region: String,
    #[arg(short, long)]
    output: String,
    #[arg(short, long)]
    annotation: Option<String>,
    #[arg(long, value_enum, default_value = "auto")]
    annotation_format: SnapAnnotationFormat,
    #[arg(long)]
    reference: Option<String>,
    #[arg(long)]
    bai: Option<String>,
    #[arg(long, default_value_t = 10)]
    mapq: u8,
    #[arg(long, default_value_t = 500)]
    max_reads: usize,
    #[arg(long, default_value_t = 1400)]
    width: u32,
    #[arg(long, default_value_t = 500)]
    min_height: u32,
    #[arg(long, default_value_t = false)]
    no_reference: bool,
    #[arg(long, default_value_t = false)]
    no_genes: bool,
    #[arg(long, default_value_t = false)]
    squash: bool,
    /// Where to draw colored base tracks: top reference strip, bottom sample strip, both, or none.
    #[arg(long, value_enum, default_value = "top")]
    base_track: SnapBaseTrackMode,
    #[arg(long, value_enum, default_value = "auto")]
    format: SnapFormat,
    /// Path to a TSV file specifying genomic markers to overlay (e.g. variants, SVs). Format: chrom\tposition\tlabel\tmarker_type\t[end_position]
    #[arg(long)]
    markers: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum RunMode {
    Scan,
    Summarize,
    Both,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum IsoformSelectMode {
    Longest,
    Shortest,
    Median,
    Common3p,
    Canonical,
}

#[derive(ClapArgs, Debug)]
struct FastqArgs {
    #[arg(short, long, num_args = 1..)]
    input: Vec<String>,
    #[arg(short, long)]
    output: Option<String>,
    #[arg(short, long, default_value_t = num_cpus::get().max(1))]
    threads: usize,
    #[arg(long, default_value_t = 1_000_000)]
    sample_size: usize,
    #[arg(long, default_value_t = 7)]
    kmer_size: usize,
    #[arg(long, default_value_t = 20)]
    top_n: usize,
    #[arg(long, default_value_t = 33)]
    phred_offset: u8,
    #[arg(long, default_value_t = false)]
    no_kmers: bool,
    #[arg(long, default_value_t = false)]
    paired: bool,
    #[arg(long, default_value_t = 1000)]
    length_bin_size: usize,
    #[arg(long, default_value_t = false)]
    pigz: bool,
    #[arg(long, default_value_t = 2)]
    pigz_threads: usize,
    #[arg(long, default_value_t = 50_000)]
    batch_size: usize,
}

#[derive(ClapArgs, Debug)]
struct AlignArgs {
    #[arg(short, long, num_args = 1..)]
    input: Vec<String>,
    #[arg(short, long)]
    output: Option<String>,
    #[arg(long, default_value_t = 10)]
    mapq: u8,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
}

#[derive(ClapArgs, Debug)]
struct RnaArgs {
    #[arg(short, long, num_args = 1..)]
    input: Vec<String>,
    #[arg(short, long, alias = "ref-bed")]
    annotation: String,
    #[arg(long, default_value = "auto")]
    annotation_format: String,
    #[arg(short, long, value_enum, default_value = "both")]
    mode: RunMode,
    #[arg(
        long,
        value_delimiter = ',',
        default_values_t = [AnalysisType::GeneBody, AnalysisType::ThreePrime, AnalysisType::Distribution, AnalysisType::Qc]
    )]
    analysis: Vec<AnalysisType>,
    #[arg(short, long)]
    output: String,
    #[arg(long, default_value_t = 200)]
    three_prime_cluster_window: u64,
    #[arg(long, default_value_t = 100)]
    min_length: u64,
    #[arg(long, default_value_t = 500)]
    normalization_bp: usize,
    #[arg(long)]
    gene_id_delimiter: Option<char>,
    #[arg(long)]
    gene_id_regex: Option<String>,
    #[arg(long, default_value_t = 10)]
    mapq: u8,
    #[arg(long, default_value_t = 15000)]
    max_3p_dist: usize,
    #[arg(long, default_value_t = 100)]
    min_support: usize,
    #[arg(long, default_value_t = false)]
    no_plot: bool,
    #[arg(long, default_value_t = false)]
    snap_qc: bool,
    #[arg(long, value_delimiter = ',')]
    snap_genes: Vec<String>,
    #[arg(long, default_value_t = 500)]
    snap_flank: u64,
    #[arg(long, default_value_t = 500)]
    snap_max_reads: usize,
    #[arg(long, default_value_t = false)]
    r2_only: bool,
    #[arg(long, default_value = "protein_coding")]
    biotype: String,
    #[arg(long, default_value = "all")]
    distribution_biotype: String,
    #[arg(long, default_value_t = false)]
    distribution_use_all_biotypes: bool,
    #[arg(long, default_value_t = false)]
    save_index: bool,
    #[arg(long, default_value_t = false)]
    load_index: bool,
    #[arg(long, default_value_t = false)]
    ends: bool,
    #[arg(long, default_value_t = false)]
    transcript_centric: bool,
    #[arg(long, default_value_t = false)]
    plus: bool,
    #[arg(long, value_enum, default_value = "common3p")]
    isoform_select: IsoformSelectMode,
    #[arg(long, default_value_t = false)]
    coverage_weighted: bool,
    #[arg(long, default_value_t = false)]
    strict_cluster: bool,
    #[arg(long, default_value_t = false)]
    stratify_length: bool,
    #[arg(long, default_value_t = 1_000_000)]
    qc_sample_size: usize,
    #[arg(long)]
    rdna_bed: Option<String>,
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "rdna,rrna,45s,18s,28s,5s,rn45s,rn18s,rn28s,rn5s"
    )]
    rdna_contigs: Vec<String>,
    #[arg(long, default_value_t = 1)]
    step_size: usize,
    #[arg(long, default_value_t = 50)]
    three_prime_bin_size: usize,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
    #[arg(long, default_value_t = false)]
    pub skip_distribution: bool,
    #[arg(long, default_value_t = 100)]
    pub three_prime_min_anchor_count: u64,
    #[arg(long, default_value_t = 5.0)]
    pub three_prime_min_anchor_mean: f64,
    #[arg(long, default_value_t = 3)]
    pub three_prime_min_anchor_nonzero_bins: usize,
    #[arg(long, default_value_t = 3.0)]
    pub three_prime_max_ratio: f64,
    #[arg(long, default_value_t = false)]
    pub write_counts: bool,
    #[arg(long, default_value_t = false)]
    pub write_gene_profiles: bool,
    #[arg(long, default_value_t = true)]
    pub compact_gene_profiles: bool,
    /// Number of concurrent threads for dense map construction (reduces peak memory).
    #[arg(long, default_value_t = 1)]
    pub dense_map_workers: usize,
    /// Strategy for building dense maps (all, chrom, chunk, or window). Window is recommended for performance and low memory.
    #[arg(long, value_enum, default_value_t = crate::analysis::rna_qc::DenseMapScope::Window)]
    pub dense_map_scope: crate::analysis::rna_qc::DenseMapScope,
    /// Size of genomic chunks for dense map construction in 'chunk' mode.
    #[arg(long, default_value_t = 50_000_000)]
    pub dense_map_chunk_size: u64,
    /// Size of genomic windows for BAM scanning. In 'window' mode, this is also the indexing granularity.
    #[arg(long, default_value_t = 10_000_000)]
    pub window_size: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    env_logger::init();

    match cli.command {
        Commands::Rna(args) => run_rna_wrapper(args),
        Commands::Fastq(args) => run_fastq(args),
        Commands::Align(args) => run_align(args),
        Commands::Dna(args) => run_dna(args),
        Commands::Atac(args) => planned_subcommand("atac", &args),
        Commands::Contam(args) => run_contam(args),
        Commands::Report(args) => run_report(args),
        Commands::Snap(args) => run_snap(args),
    }
}

fn run_rna_wrapper(args: RnaArgs) -> Result<()> {
    let mut analysis = args.analysis.clone();
    if args.skip_distribution {
        analysis.retain(|a| *a != AnalysisType::Distribution);
    }

    let config = RnaQcConfig {
        input: args.input,
        annotation: args.annotation,
        annotation_format: args.annotation_format.clone(),
        output: args.output,
        threads: args.threads,
        analysis,
        three_prime_cluster_window: args.three_prime_cluster_window,
        min_length: args.min_length,
        normalization_bp: args.normalization_bp,
        gene_id_delimiter: args.gene_id_delimiter,
        gene_id_regex: args.gene_id_regex,
        mapq: args.mapq,
        max_3p_dist: args.max_3p_dist,
        min_support: args.min_support,
        no_plot: args.no_plot,
        snap_qc: args.snap_qc,
        snap_genes: args.snap_genes,
        snap_flank: args.snap_flank,
        snap_max_reads: args.snap_max_reads,
        r2_only: args.r2_only,
        biotype: args.biotype,
        distribution_biotype: args.distribution_biotype,
        distribution_use_all_biotypes: args.distribution_use_all_biotypes,
        save_index: args.save_index,
        load_index: args.load_index,
        ends: args.ends,
        transcript_centric: args.transcript_centric,
        plus: args.plus,
        isoform_select: match args.isoform_select {
            IsoformSelectMode::Longest => crate::io::annotation::IsoformSelect::Longest,
            IsoformSelectMode::Shortest => crate::io::annotation::IsoformSelect::Shortest,
            IsoformSelectMode::Median => crate::io::annotation::IsoformSelect::Median,
            IsoformSelectMode::Common3p => crate::io::annotation::IsoformSelect::Common3p,
            IsoformSelectMode::Canonical => crate::io::annotation::IsoformSelect::Canonical,
        },
        coverage_weighted: args.coverage_weighted,
        strict_cluster: args.strict_cluster,
        stratify_length: args.stratify_length,
        qc_sample_size: args.qc_sample_size,
        rdna_bed: args.rdna_bed,
        rdna_contigs: args.rdna_contigs,
        step_size: args.step_size,
        three_prime_bin_size: args.three_prime_bin_size,
        three_prime_min_anchor_count: args.three_prime_min_anchor_count,
        three_prime_min_anchor_mean: args.three_prime_min_anchor_mean,
        three_prime_min_anchor_nonzero_bins: args.three_prime_min_anchor_nonzero_bins,
        three_prime_max_ratio: args.three_prime_max_ratio,
        write_counts: args.write_counts,
        write_gene_profiles: args.write_gene_profiles,
        compact_gene_profiles: args.compact_gene_profiles,
        dense_map_workers: args.dense_map_workers,
        dense_map_scope: args.dense_map_scope,
        dense_map_chunk_size: args.dense_map_chunk_size,
        window_size: args.window_size,
    };
    run_rna(config)
}

fn run_align(args: AlignArgs) -> Result<()> {
    let output_prefix = args.output.unwrap_or_else(|| {
        if args.input.len() == 1 {
            sample_name_from_alignment_path(&args.input[0])
        } else {
            "rs-qc".to_string()
        }
    });

    println!("--------------------------------------------------");
    println!("rs-qc align v{}", env!("CARGO_PKG_VERSION"));
    println!("--------------------------------------------------");
    println!("Scanning {} alignment input(s).", args.input.len());

    let config = AlignmentQcConfig {
        inputs: args.input,
        output_prefix: output_prefix.clone(),
        mapq_threshold: args.mapq,
        threads: args.threads,
        show_progress: true,
    };

    let metrics = run_alignment_qc(&config)?;
    println!("  - Records: {}", metrics.total_records);
    println!(
        "  - Mapped: {} ({:.2}%)",
        metrics.mapped_records,
        percent(metrics.mapped_records, metrics.total_records)
    );
    println!(
        "  - Properly paired: {} ({:.2}% of paired records)",
        metrics.properly_paired_records,
        percent(metrics.properly_paired_records, metrics.paired_records)
    );
    println!(
        "  - de:f accuracy records: {}, median: {}",
        metrics.accuracy.records_with_de,
        metrics
            .accuracy
            .median()
            .map(|v| format!("{v:.4}"))
            .unwrap_or_else(|| "NA".to_string())
    );
    if let Some(ref m) = metrics.insert_size_metrics {
        println!("  - Insert Size & cfDNA Diagnostics:");
        println!(
            "    * Mono-nucleosomal Peak: {} bp",
            m.mono_nucleosomal_peak.map(|x| x.to_string()).unwrap_or_else(|| "NA".to_string())
        );
        println!(
            "    * Di-nucleosomal Peak:   {} bp",
            m.di_nucleosomal_peak.map(|x| x.to_string()).unwrap_or_else(|| "NA".to_string())
        );
        println!(
            "    * cfDNA Short/Mono Ratio: {}",
            m.cfdna_ratio.map(|v| format!("{v:.4}")).unwrap_or_else(|| "NA".to_string())
        );
        println!(
            "    * Mono/Di Ratio:          {}",
            m.mono_di_ratio.map(|v| format!("{v:.4}")).unwrap_or_else(|| "NA".to_string())
        );
        println!("    * Sub-nucleosomal Frac:   {:.2}%", m.short_fraction * 100.0);
    }
    println!(
        "  - Summary written to: {}.align.summary.txt",
        output_prefix
    );
    Ok(())
}

fn run_fastq(args: FastqArgs) -> Result<()> {
    let output_prefix = args.output.unwrap_or_else(|| {
        if args.input.len() == 1 {
            sample_name_from_path(&args.input[0])
        } else {
            "rs-qc".to_string()
        }
    });

    println!("--------------------------------------------------");
    println!("rs-qc fastq v{}", env!("CARGO_PKG_VERSION"));
    println!("--------------------------------------------------");
    println!("Scanning {} FASTQ input(s).", args.input.len());

    let config = FastqQcConfig {
        inputs: args.input,
        output_prefix: output_prefix.clone(),
        threads: args.threads.max(1),
        sample_size: args.sample_size,
        kmer_size: args.kmer_size,
        top_n: args.top_n,
        phred_offset: args.phred_offset,
        no_kmers: args.no_kmers,
        paired: args.paired,
        length_bin_size: args.length_bin_size,
        use_pigz: args.pigz,
        pigz_threads: args.pigz_threads.max(1),
        batch_size: args.batch_size.max(1),
    };

    let metrics = run_fastq_qc(&config)?;
    println!("  - Reads: {}", metrics.total_reads);
    println!("  - Bases: {}", metrics.total_bases);
    println!("  - Mean read length: {:.2}", metrics.mean_read_length());
    if let Some(n50) = metrics.read_nx().get(&50) {
        println!("  - Read N50: {}", n50);
    }
    println!(
        "  - Summary written to: {}.fastq.summary.txt / .json",
        output_prefix
    );
    Ok(())
}

fn run_dna(args: DnaArgs) -> Result<()> {
    let output_prefix = args.output.clone().unwrap_or_else(|| {
        if args.input.len() == 1 {
            sample_name_from_alignment_path(&args.input[0])
        } else {
            "rs-qc".to_string()
        }
    });

    println!("--------------------------------------------------");
    println!("rs-qc dna v{}", env!("CARGO_PKG_VERSION"));
    println!("--------------------------------------------------");
    println!("Scanning {} BAM input(s).", args.input.len());

    let config = DnaQcConfig {
        inputs: args.input,
        output_prefix: output_prefix.clone(),
        mapq_threshold: args.mapq,
        threads: args.threads,
        window_size: args.window_size,
        targets_path: args.targets,
        thresholds: args.thresholds,
        callable_depth: args.callable_depth,
        include_duplicates: args.include_duplicates,
        reference_fasta: args.reference,
        annotation_path: args.annotation,
        low_cov_threshold: args.low_cov_threshold,
        snap_lowcov: args.snap_lowcov,
        show_progress: true,
    };

    run_dna_qc(&config)?;
    println!(
        "  - Summary written to: {}.dna.summary.txt / .json",
        output_prefix
    );
    Ok(())
}

fn run_contam(args: ContamArgs) -> Result<()> {
    let output_prefix = args.output.unwrap_or_else(|| {
        if args.input.len() == 1 {
            sample_name_from_alignment_path(&args.input[0])
        } else {
            "rs-qc".to_string()
        }
    });

    println!("--------------------------------------------------");
    println!("rs-qc contam v{}", env!("CARGO_PKG_VERSION"));
    println!("--------------------------------------------------");
    println!("Scanning {} BAM input(s).", args.input.len());
    if let Some(path) = &args.rrna_fasta {
        println!(
            "  - Exact rRNA k-mer screen: {} (k={}, min hits={})",
            path, args.kmer_size, args.min_kmer_hits
        );
    }

    let config = ContaminationQcConfig {
        inputs: args.input,
        output_prefix: output_prefix.clone(),
        rdna_bed: args.rdna_bed,
        rdna_contigs: args.rdna_contigs,
        rrna_fasta: args.rrna_fasta,
        kmer_size: args.kmer_size,
        min_kmer_hits: args.min_kmer_hits,
        sample_size: args.sample_size,
        mapq_threshold: args.mapq,
        scan_all_reads: args.kmer_scan_all_reads,
    };

    let metrics = run_contamination_qc(&config)?;
    for summary in metrics {
        println!("\nSample: {}", summary.sample);
        println!(
            "    rDNA: interval={} reads, contig={} reads ({:.2}%)",
            summary.rdna_interval_reads,
            summary.rdna_contig_reads,
            summary.rdna_aligned_fraction * 100.0
        );
        if summary.rrna_kmer_reads_screened > 0 {
            println!(
                "    rRNA k-mer: {} / {} screened ({:.2}%)",
                summary.rrna_kmer_reads,
                summary.rrna_kmer_reads_screened,
                summary.rrna_kmer_fraction * 100.0
            );
        }
    }
    println!(
        "  - Summary written to: {}.<sample>.contam.summary.tsv / .json",
        output_prefix
    );
    Ok(())
}

fn run_report(args: ReportArgs) -> Result<()> {
    let output_prefix = args.output.unwrap_or_else(|| "rs-qc-report".to_string());
    println!("--------------------------------------------------");
    println!("rs-qc report v{}", env!("CARGO_PKG_VERSION"));
    println!("--------------------------------------------------");
    println!("Collecting summary JSON files...");

    let sources = resolve_input_files(&args.input)?;
    println!("  - Found {} summary file(s).", sources.len());

    let document = build_document(&sources)?;
    write_report(&document, &output_prefix)?;
    println!("  - Report written to: {}.html", output_prefix);
    Ok(())
}

fn run_snap(args: SnapArgs) -> Result<()> {
    let region = resolve_snapshot_region(
        &args.region,
        args.annotation.as_deref(),
        match args.annotation_format {
            SnapAnnotationFormat::Auto => AnnotationFormat::Auto,
            SnapAnnotationFormat::Gtf => AnnotationFormat::Gtf,
            SnapAnnotationFormat::Bed12 => AnnotationFormat::Bed12,
        },
    )?;
    validate_output_format(
        &args.output,
        match args.format {
            SnapFormat::Auto => SnapshotOutputFormat::Auto,
            SnapFormat::Png => SnapshotOutputFormat::Png,
            SnapFormat::Svg => SnapshotOutputFormat::Svg,
        },
    )?;

    println!("--------------------------------------------------");
    println!("rs-qc snap v{}", env!("CARGO_PKG_VERSION"));
    println!("--------------------------------------------------");
    println!(
        "Rendering {}:{}-{}",
        region.chrom,
        region.start + 1,
        region.end
    );

    let (show_reference_base_track, show_sample_base_track) = match args.base_track {
        SnapBaseTrackMode::Top => (true, false),
        SnapBaseTrackMode::Bottom => (false, true),
        SnapBaseTrackMode::Both => (true, true),
        SnapBaseTrackMode::None => (false, false),
    };

    let config = SnapshotConfig {
        bam_path: args.bam,
        bai_path: args.bai,
        region,
        annotation_path: args.annotation,
        annotation_format: match args.annotation_format {
            SnapAnnotationFormat::Auto => AnnotationFormat::Auto,
            SnapAnnotationFormat::Gtf => AnnotationFormat::Gtf,
            SnapAnnotationFormat::Bed12 => AnnotationFormat::Bed12,
        },
        reference_path: args.reference,
        output_path: args.output.clone(),
        mapq_threshold: args.mapq,
        max_reads: args.max_reads,
        width: args.width,
        min_height: args.min_height,
        show_reference: !args.no_reference,
        show_genes: !args.no_genes,
        show_reference_base_track,
        show_sample_base_track,
        squash: args.squash,
        markers_path: args.markers,
        inline_markers: Vec::new(),
    };

    run_snapshot(&config)?;
    println!("  - Snapshot written to: {}", args.output);
    Ok(())
}

fn percent(n: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        n as f64 * 100.0 / total as f64
    }
}

fn planned_subcommand(name: &str, args: &PlannedArgs) -> Result<()> {
    bail!(
        "rs-qc {name} is planned but not implemented yet (output prefix: {}). Use `rs-qc rna ...` for the current RNA-seq QC pipeline.",
        args.output
    )
}
