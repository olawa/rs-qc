mod analysis;
mod io;
mod models;
mod stats;

use crate::analysis::alignment_qc::{
    run_alignment_qc, sample_name_from_alignment_path, AlignmentQcConfig,
};
use crate::analysis::bam_scan::{find_bai_path, generate_windows};
use crate::analysis::contamination::ContaminantIndex;
use crate::analysis::fastq_qc::{run_fastq_qc, sample_name_from_path, FastqQcConfig};
use crate::analysis::feature_index::FeatureIndex;
use crate::analysis::index::{AnnotationIndex, DenseMap};
use crate::analysis::qc::{InlineQcState, ReadEndObservation};
use crate::analysis::read_distribution::RegionType;
use crate::analysis::types::AnalysisType;
use crate::io::annotation::{load_annotation, load_genes, AnnotationConfig, AnnotationFormat, IsoformSelect};
use crate::models::{normalize_chrom, Gene};
use crate::stats::plotting::{
    generate_gene_body_plot, generate_multi_3p_dist_plot, generate_multi_inner_distance_plot,
    generate_stratified_gene_body_plot,
};
use crate::stats::{
    aggregate_genes, aggregate_rseqc_classic, aggregate_rseqc_stratified,
    write_classic_wide_format,
};
use anyhow::{bail, Result};
use clap::{Args as ClapArgs, Parser, Subcommand, ValueEnum};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use noodles::bam;
// use rayon::prelude::*; // redundant because of local use in main loop
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::sync::Arc;
use std::time::Instant;

thread_local! {
    static BAM_READER: RefCell<Option<bam::io::IndexedReader<noodles::bgzf::Reader<File>>>> = RefCell::new(None);
}

const MAX_REASONABLE_INNER_DISTANCE_BP: i64 = 2_000;

struct LocalGeneBuffers {
    diff_3p: Vec<i32>,
    diff_percentile: Vec<i32>,
}

impl LocalGeneBuffers {
    fn new(max_3p_dist: usize) -> Self {
        Self {
            diff_3p: vec![0; max_3p_dist + 1],
            diff_percentile: vec![0; 101],
        }
    }

    fn add(&mut self, s_5p_idx: usize, total_len: usize) {
        let pct_idx = (s_5p_idx * 100) / total_len.max(1);
        if pct_idx < 100 {
            self.diff_percentile[pct_idx] += 1;
            self.diff_percentile[pct_idx + 1] -= 1;
        }

        let dist_3p = total_len.saturating_sub(s_5p_idx + 1);
        if dist_3p < self.diff_3p.len() - 1 {
            self.diff_3p[dist_3p] += 1;
            self.diff_3p[dist_3p + 1] -= 1;
        }
    }
}

struct WorkerState {
    records_seen: u64,
    fail_unmapped: u64,
    fail_secondary: u64,
    fail_qc: u64,
    fail_mapq: u64,
    overlaps_found: u64,
    total_tags: u64,
    aligned_qc_reads: u64,
    mtdna_reads: u64,
    rdna_reads: u64,
    read_dist_counts: HashMap<RegionType, u64>,
    /// Flat array indexed by gene_idx. None = gene not yet seen this worker.
    /// Avoids HashMap overhead (hash + probe) on every base in the hot loop.
    gene_counts: Vec<Option<LocalGeneBuffers>>,
    /// Gene indices whose slot in gene_counts is Some — used to restrict merge
    /// to only live entries without scanning the full Vec.
    populated: Vec<u32>,
    max_3p_dist: usize,
    qc: InlineQcState,
}

impl WorkerState {
    fn new_with_capacity(n_genes: usize, max_3p_dist: usize, qc_sample_size: usize) -> Self {
        Self {
            records_seen: 0,
            fail_unmapped: 0,
            fail_secondary: 0,
            fail_qc: 0,
            fail_mapq: 0,
            overlaps_found: 0,
            total_tags: 0,
            aligned_qc_reads: 0,
            mtdna_reads: 0,
            rdna_reads: 0,
            read_dist_counts: HashMap::new(),
            gene_counts: (0..n_genes).map(|_| None).collect(),
            populated: Vec::new(),
            max_3p_dist,
            qc: InlineQcState::new(qc_sample_size),
        }
    }

    /// Return a mutable reference to the buffer for `g_idx`, lazily initialising
    /// it if this is the first time this worker has seen the gene.
    #[inline(always)]
    fn get_or_init(&mut self, g_idx: usize) -> &mut LocalGeneBuffers {
        if self.gene_counts[g_idx].is_none() {
            self.gene_counts[g_idx] = Some(LocalGeneBuffers::new(self.max_3p_dist));
            self.populated.push(g_idx as u32);
        }
        self.gene_counts[g_idx].as_mut().unwrap()
    }

    fn merge(mut self, other: Self) -> Self {
        self.records_seen += other.records_seen;
        self.fail_unmapped += other.fail_unmapped;
        self.fail_secondary += other.fail_secondary;
        self.fail_qc += other.fail_qc;
        self.fail_mapq += other.fail_mapq;
        self.overlaps_found += other.overlaps_found;
        self.total_tags += other.total_tags;
        self.aligned_qc_reads += other.aligned_qc_reads;
        self.mtdna_reads += other.mtdna_reads;
        self.rdna_reads += other.rdna_reads;

        for (region, count) in other.read_dist_counts {
            *self.read_dist_counts.entry(region).or_insert(0) += count;
        }

        self.qc.merge_from(other.qc);

        // Only touch genes that other actually wrote to.
        for g_idx in other.populated {
            let g = g_idx as usize;
            if let Some(other_buf) = &other.gene_counts[g] {
                let entry = self.get_or_init(g);
                for (i, &val) in other_buf.diff_3p.iter().enumerate() {
                    entry.diff_3p[i] += val;
                }
                for (i, &val) in other_buf.diff_percentile.iter().enumerate() {
                    entry.diff_percentile[i] += val;
                }
            }
        }
        self
    }
}

#[derive(ValueEnum, Clone, Debug, PartialEq, Eq, Default)]
pub enum RunMode {
    RseqcClassic,
    ThreePrime,
    #[default]
    Both,
}

/// CLI-facing mirror of `IsoformSelect` — kept separate to avoid adding clap
/// as a dependency of the annotation module.
#[derive(ValueEnum, Clone, Debug, PartialEq, Eq, Default)]
pub enum IsoformSelectMode {
    /// Longest spliced transcript (classic RSeQC behaviour).
    #[default]
    Longest,
    /// Shortest qualifying transcript. Reduces 5' bias caused by long un-expressed 3' UTRs.
    Shortest,
    /// Transcript closest to the median length — a compromise between the two.
    Median,
}

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
    Dna(PlannedArgs),
    /// ATAC/ChIP/cfDNA-style assay metrics.
    Atac(PlannedArgs),
    /// Contamination metrics from aligned intervals, contigs, and optional k-mer screens.
    Contam(PlannedArgs),
    /// Render JSON/TSV/SVG metrics into a unified report.
    Report(PlannedArgs),
}

#[derive(ClapArgs, Debug)]
struct PlannedArgs {
    /// Output prefix for future module outputs.
    #[arg(short, long, default_value = "rs-qc")]
    output: String,
}

#[derive(ClapArgs, Debug)]
struct FastqArgs {
    #[arg(short, long, num_args = 1..)]
    input: Vec<String>,
    #[arg(short, long)]
    output: Option<String>,
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
    #[arg(long, default_value_t = 10000)]
    max_3p_dist: usize,
    #[arg(long, default_value_t = 100)]
    min_support: usize,
    #[arg(long, default_value_t = false)]
    no_plot: bool,
    #[arg(long, default_value_t = false)]
    r2_only: bool,
    #[arg(long, default_value = "protein_coding")]
    biotype: String,
    #[arg(long, default_value_t = false)]
    save_index: bool,
    #[arg(long, default_value_t = false)]
    load_index: bool,
    #[arg(long, default_value_t = false)]
    ends: bool,
    #[arg(long, default_value_t = false)]
    transcript_centric: bool,
    /// Only load transcripts on the '+' strand. Halves the index size and avoids
    /// mixing orientation artifacts. Recommended for unstranded or FR-stranded libraries.
    #[arg(long, default_value_t = false)]
    plus: bool,
    /// Which isoform to use as the per-gene representative in the 100-bin gene body
    /// coverage plot. 'longest' matches classic RSeQC behaviour but can introduce 5'
    /// bias when long 3' UTRs are unannotated or not expressed. 'shortest' or 'median'
    /// produce profiles closer to curated BED-file results.
    #[arg(long, value_enum, default_value = "longest")]
    isoform_select: IsoformSelectMode,
    /// Weight each gene's contribution to the 100-bin plot by its coverage depth
    /// (capped at 20× min_support). Off by default to match RSeQC behaviour where
    /// every expressed gene contributes equally after per-gene normalisation.
    #[arg(long, default_value_t = false)]
    coverage_weighted: bool,
    /// Exclude genes where the spread of isoform 3' ends exceeds
    /// `--three-prime-cluster-window`. Set --three-prime-cluster-window (default 200 bp)
    /// to your desired maximum spread, e.g. 1000 for 1 kb. Genes with highly divergent
    /// alternative polyadenylation sites are removed, reducing 5' artefacts in the
    /// gene body coverage plot.
    #[arg(long, default_value_t = false)]
    strict_cluster: bool,
    /// Additionally produce a stratified gene body coverage SVG where genes are split
    /// by transcript length: short (<1.5 kb), medium (1.5–5 kb), long (>5 kb).
    /// Useful for diagnosing whether coverage artefacts are length-dependent.
    #[arg(long, default_value_t = false)]
    stratify_length: bool,
    #[arg(long, default_value_t = 100_000)]
    qc_sample_size: usize,
    #[arg(long)]
    rdna_bed: Option<String>,
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "rdna,rrna,45s,18s,28s,5s,rn45s,rn18s,rn28s,rn5s"
    )]
    rdna_contigs: Vec<String>,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
}

fn normalized_name_set(names: &[String]) -> HashSet<String> {
    names.iter().map(|name| normalize_chrom(name)).collect()
}

fn is_mtdna_chrom(chrom: &str) -> bool {
    matches!(normalize_chrom(chrom).as_str(), "m" | "mt" | "mitochondria")
}

fn is_rdna_chrom(chrom: &str, rdna_contigs: &HashSet<String>) -> bool {
    let chrom = normalize_chrom(chrom);
    rdna_contigs.contains(&chrom)
        || chrom.contains("rdna")
        || chrom.contains("rrna")
        || chrom.contains("ribosomal")
}

fn aligned_match_span(record: &bam::Record) -> Option<(u64, u64)> {
    let start = match record.alignment_start() {
        Some(Ok(p)) => p.get() as u64 - 1,
        _ => return None,
    };

    let mut curr = start;
    let mut first = None;
    let mut last = None;

    for op_res in record.cigar().iter() {
        if let Ok(op) = op_res {
            use noodles::sam::alignment::record::cigar::op::Kind;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if first.is_none() {
                        first = Some(curr);
                    }
                    last = Some(curr + op.len() as u64 - 1);
                    curr += op.len() as u64;
                }
                Kind::Deletion | Kind::Skip => {
                    curr += op.len() as u64;
                }
                _ => {}
            }
        }
    }

    first.zip(last)
}

fn maybe_observe_inline_qc(
    state: &mut InlineQcState,
    record: &bam::Record,
    chrom: &str,
    dense: &crate::analysis::index::DenseMap,
    genes: &[Gene],
    win_start: u64,
    win_end: u64,
) {
    let flags = record.flags();
    if !state.needs_more() || !flags.is_segmented() {
        return;
    }
    if !flags.is_first_segment() && !flags.is_last_segment() {
        return;
    }

    let aln_start = match record.alignment_start() {
        Some(Ok(p)) => p.get() as u64 - 1,
        _ => return,
    };
    if aln_start < win_start || aln_start >= win_end {
        return;
    }

    if matches!(
        dense.get_hits(aln_start),
        crate::analysis::index::Hits::None
    ) {
        return;
    }

    let (first_match, last_match) = match aligned_match_span(record) {
        Some(span) => span,
        None => return,
    };

    let qname = match record.name() {
        Some(name) => name.as_ref().to_vec(),
        None => return,
    };

    let observed = ReadEndObservation {
        chrom: chrom.to_string(),
        first_match,
        last_match,
        is_first: flags.is_first_segment(),
        is_reverse: flags.is_reverse_complemented(),
    };

    let Some((read1, read2)) = state.observe_end(qname, observed) else {
        return;
    };

    if read1.chrom != read2.chrom || read1.chrom != chrom {
        return;
    }

    let gene_idx = match (
        dense.get_hits(read1.first_match),
        dense.get_hits(read2.first_match),
    ) {
        (
            crate::analysis::index::Hits::Single(left_idx),
            crate::analysis::index::Hits::Single(right_idx),
        ) if left_idx == right_idx => left_idx as usize,
        _ => return,
    };

    let gene = &genes[gene_idx];
    if gene.representative.strand != '+' {
        return;
    }

    let tx_start_1 = match gene.bin_map.get_spliced_5p(read1.first_match) {
        Some(pos) => pos,
        None => return,
    };
    let tx_end_1 = match gene.bin_map.get_spliced_5p(read1.last_match) {
        Some(pos) => pos + 1,
        None => return,
    };
    let tx_start_2 = match gene.bin_map.get_spliced_5p(read2.first_match) {
        Some(pos) => pos,
        None => return,
    };
    let tx_end_2 = match gene.bin_map.get_spliced_5p(read2.last_match) {
        Some(pos) => pos + 1,
        None => return,
    };

    let inner_distance = tx_start_1.max(tx_start_2) as i64 - tx_end_1.min(tx_end_2) as i64;
    if inner_distance.abs() > MAX_REASONABLE_INNER_DISTANCE_BP {
        return;
    }

    state.record_pair(inner_distance, read1.is_reverse, read2.is_reverse);
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    env_logger::init();

    match cli.command {
        Commands::Rna(args) => run_rna(args),
        Commands::Fastq(args) => run_fastq(args),
        Commands::Align(args) => run_align(args),
        Commands::Dna(args) => planned_subcommand("dna", &args),
        Commands::Atac(args) => planned_subcommand("atac", &args),
        Commands::Contam(args) => planned_subcommand("contam", &args),
        Commands::Report(args) => planned_subcommand("report", &args),
    }
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
        sample_size: args.sample_size,
        kmer_size: args.kmer_size,
        top_n: args.top_n,
        phred_offset: args.phred_offset,
        no_kmers: args.no_kmers,
        paired: args.paired,
    };

    let metrics = run_fastq_qc(&config)?;
    println!("  - Reads: {}", metrics.total_reads);
    println!("  - Bases: {}", metrics.total_bases);
    println!("  - Mean read length: {:.2}", metrics.mean_read_length());
    println!(
        "  - Sampled duplication estimate: {:.2}%",
        metrics.duplication_estimate() * 100.0
    );
    if config.paired {
        println!(
            "  - Paired reads checked: {}, name mismatches: {}",
            metrics.paired_reads_checked, metrics.paired_name_mismatches
        );
    }
    println!(
        "  - Summary written to: {}.fastq.summary.txt",
        output_prefix
    );
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

fn run_rna(args: RnaArgs) -> Result<()> {
    let threads = args.threads;
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let start_time = Instant::now();

    let ann_config = AnnotationConfig {
        gene_id_delimiter: args.gene_id_delimiter,
        gene_id_regex: args.gene_id_regex,
        biotype_filter: Some(args.biotype.clone()),
        three_prime_cluster_window: args.three_prime_cluster_window,
        min_transcript_length: args.min_length,
        max_3p_dist: args.max_3p_dist,
        transcript_centric: args.transcript_centric,
        plus_strand_only: args.plus,
        isoform_select: match args.isoform_select {
            IsoformSelectMode::Longest  => IsoformSelect::Longest,
            IsoformSelectMode::Shortest => IsoformSelect::Shortest,
            IsoformSelectMode::Median   => IsoformSelect::Median,
        },
        strict_cluster: args.strict_cluster,
    };
    let ann_format = match args.annotation_format.to_lowercase().as_str() {
        "bed12" => AnnotationFormat::Bed12,
        "gtf" => AnnotationFormat::Gtf,
        _ => AnnotationFormat::Auto,
    };

    println!("--------------------------------------------------");
    println!("rs-qc rna v{}", env!("CARGO_PKG_VERSION"));
    println!("--------------------------------------------------");
    println!("Loading annotation: {}", args.annotation);

    let index_path = format!("{}.ridx", args.annotation);
    let dist_enabled = args.analysis.contains(&AnalysisType::Distribution);
    let (index_proto, feature_index) = if dist_enabled {
        if args.load_index {
            println!(
                "  - Distribution analysis needs transcript annotations, so the source annotation will be parsed instead of loading only the collapsed index."
            );
        }
        let loaded = load_annotation(&args.annotation, ann_format, &ann_config)?;
        println!("  - Extracted {} eligible genes.", loaded.genes.len());
        let index = AnnotationIndex::new(loaded.genes, false);
        let feature_index = Arc::new(FeatureIndex::build(&loaded.transcripts));
        if args.save_index {
            println!("  - Saving index for future use: {}", index_path);
            index.save_to_file(&index_path)?;
        }
        (index, Some(feature_index))
    } else if args.load_index && std::path::Path::new(&index_path).exists() {
        println!("  - Loading pre-built index: {}", index_path);
        (
            AnnotationIndex::load_from_file(&index_path, args.max_3p_dist)?,
            None,
        )
    } else {
        let genes = load_genes(&args.annotation, ann_format, &ann_config)?;
        println!("  - Extracted {} eligible genes.", genes.len());
        let index = AnnotationIndex::new(genes, false);

        if args.save_index {
            println!("  - Saving index for future use: {}", index_path);
            index.save_to_file(&index_path)?;
        }
        (index, None)
    };

    println!(
        "  - Index ready with {} active genes.",
        index_proto.genes.len()
    );

    let rdna_contigs = Arc::new(normalized_name_set(&args.rdna_contigs));
    let rdna_intervals = if let Some(path) = &args.rdna_bed {
        println!("  - Loading rDNA intervals: {}", path);
        Some(Arc::new(ContaminantIndex::from_bed(path)?))
    } else {
        None
    };

    let mut classic_all = HashMap::new();
    let mut classic_percent_all = HashMap::new();
    let mut dist_3p_all = HashMap::new();
    let mut inner_distance_all = HashMap::new();
    let mut percentile_all = HashMap::new();

    for (i, bam_path) in args.input.iter().enumerate() {
        let sample_name = std::path::Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(bam_path);
        println!(
            "\n[{}/{}] Processing BAM: {}",
            i + 1,
            args.input.len(),
            bam_path
        );

        // Build dense maps for all chromosomes in this BAM
        let header = {
            let file = File::open(bam_path)?;
            let mut reader = bam::io::Reader::new(file);
            reader.read_header()?
        };

        println!("  - Pre-calculating dense mappings for active chromosomes...");
        use rayon::prelude::*;
        let dense_maps: HashMap<String, Arc<DenseMap>> = header
            .reference_sequences()
            .iter()
            .collect::<Vec<_>>()
            .into_par_iter()
            .filter_map(|(name, _target)| {
                let chrom = String::from_utf8_lossy(name.as_ref()).to_string();
                let needs_dense = args.analysis.contains(&AnalysisType::GeneBody)
                    || args.analysis.contains(&AnalysisType::ThreePrime)
                    || args.analysis.contains(&AnalysisType::Qc);

                let dense = if needs_dense {
                    index_proto.build_dense_map(&chrom).map(Arc::new)
                } else {
                    None
                };

                if dense.is_some() {
                    Some((chrom, dense))
                } else {
                    None
                }
            })
            .fold(
                || HashMap::new(),
                |mut d_acc, (chrom, d)| {
                    if let Some(dense) = d {
                        d_acc.insert(chrom, dense);
                    }
                    d_acc
                },
            )
            .reduce(
                || HashMap::new(),
                |mut d1, d2| {
                    d1.extend(d2);
                    d1
                },
            );

        println!("  - Built maps for {} chromosomes.", dense_maps.len());
        let dense_maps = Arc::new(dense_maps);

        let mb_index_path = find_bai_path(bam_path);

        if let Some(bai_path) = mb_index_path {
            println!("  - BAI Index found. Using high-performance parallel dense scan.");
            let bai = bam::bai::read(&bai_path)?;
            let windows = generate_windows(&header, 5_000_000); // Smaller windows for better parallelization

            let m_pb = MultiProgress::new();
            let pb = m_pb.add(ProgressBar::new(windows.len() as u64));
            pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Windows ({eta})")?);

            let bai_arc = Arc::new(bai);
            let n_genes = index_proto.genes.len();
            let max_3p_dist = args.max_3p_dist;
            let qc_sample_size = args.qc_sample_size.div_ceil(threads);
            let total_state = windows
                .par_iter()
                .fold(
                    move || WorkerState::new_with_capacity(n_genes, max_3p_dist, qc_sample_size),
                    |mut state, window| {
                    BAM_READER.with(|cell| {
                        let mut opt = cell.borrow_mut();
                        if opt.is_none() {
                            let f = File::open(bam_path).expect("Failed to open BAM");
                            let mut r = bam::io::indexed_reader::Builder::default()
                                .set_index(bai_arc.as_ref().clone())
                                .build_from_reader(f)
                                .expect("Failed to build indexed reader");
                            let _ = r.read_header().expect("Failed to read BAM header");
                            *opt = Some(r);
                        }

                        let reader = opt.as_mut().unwrap();
                        let chrom = &window.chrom;
                        let maybe_dense = dense_maps.get(chrom);
                        let maybe_feature = feature_index
                            .as_ref()
                            .and_then(|index| index.chroms.get(chrom));
                        let is_mtdna_window = is_mtdna_chrom(chrom);
                        let is_rdna_contig_window = is_rdna_chrom(chrom, rdna_contigs.as_ref());
                        let maybe_rdna_intervals = rdna_intervals
                            .as_ref()
                            .and_then(|index| index.chroms.get(&normalize_chrom(chrom)));

                        let region: noodles::core::Region =
                            match format!("{}:{}-{}", chrom, window.start + 1, window.end).parse() {
                                Ok(r) => r,
                                Err(_) => return, // Can't return state from closure here easily
                            };

                        let win_start = window.start as u64;
                        let win_end = window.end as u64;
                        let mut feature_cursor = maybe_feature.map(|chrom_index| chrom_index.cursor_at(win_start));
                        let mut rdna_cursor = maybe_rdna_intervals
                            .map(|chrom_index| chrom_index.cursor_at(win_start));

                        match reader.query(&header, &region) {
                            Ok(query) => {
                                for result in query {
                                    state.records_seen += 1;
                                    match result {
                                        Ok(record) => {
                                            let flags = record.flags();

                                            if flags.is_unmapped() {
                                                state.fail_unmapped += 1;
                                                continue;
                                            }
                                            if flags.is_secondary() || flags.is_supplementary() {
                                                state.fail_secondary += 1;
                                                continue;
                                            }
                                            if flags.is_qc_fail() || flags.is_duplicate() {
                                                state.fail_qc += 1;
                                                continue;
                                            }
                                            let record_mapq = record
                                                .mapping_quality()
                                                .map(|m| m.get())
                                                .unwrap_or(255);
                                            if record_mapq < args.mapq {
                                                state.fail_mapq += 1;
                                                continue;
                                            }

                                            if args.r2_only && !flags.is_last_segment() {
                                                continue;
                                            }

                                            let pos = match record.alignment_start() {
                                                Some(Ok(p)) => p.get() as u32 - 1,
                                                _ => continue,
                                            };

                                            let genes_ref = &index_proto.genes;

                                            if (pos as u64) >= win_start && (pos as u64) < win_end {
                                                state.aligned_qc_reads += 1;
                                                if is_mtdna_window {
                                                    state.mtdna_reads += 1;
                                                }

                                                let mut is_rdna = is_rdna_contig_window;
                                                if !is_rdna {
                                                    if let Some(chrom_index) = maybe_rdna_intervals {
                                                        if let Some(cursor) = rdna_cursor.as_mut() {
                                                            if let Some((first_match, last_match)) =
                                                                aligned_match_span(&record)
                                                            {
                                                                let mid = first_match
                                                                    + (last_match - first_match) / 2;
                                                                is_rdna = chrom_index.contains(mid, cursor);
                                                            }
                                                        }
                                                    }
                                                }
                                                if is_rdna {
                                                    state.rdna_reads += 1;
                                                }
                                            }

                                            // --- Read Distribution logic ---
                                            // Classify each read once using the midpoint of its aligned span.
                                            if let Some(chrom_index) = maybe_feature {
                                                if let Some((first_match, last_match)) =
                                                    aligned_match_span(&record)
                                                {
                                                    let mid = first_match + (last_match - first_match) / 2;
                                                    if mid >= win_start && mid < win_end {
                                                        state.total_tags += 1;
                                                        if let Some(cursor) = feature_cursor.as_mut() {
                                                            let region = chrom_index.classify(mid, cursor);
                                                            *state
                                                                .read_dist_counts
                                                                .entry(region)
                                                                .or_insert(0) += 1;
                                                        }
                                                    }
                                                }
                                            }

                                            // --- Gene Body Coverage & QC logic ---
                                            if let Some(dense) = maybe_dense {
                                                if args.analysis.contains(&AnalysisType::Qc) {
                                                    maybe_observe_inline_qc(
                                                        &mut state.qc,
                                                        &record,
                                                        chrom,
                                                        dense,
                                                        genes_ref,
                                                        win_start,
                                                        win_end,
                                                    );
                                                }

                                                let needs_coverage =
                                                    args.analysis.contains(&AnalysisType::GeneBody)
                                                        || args
                                                            .analysis
                                                            .contains(&AnalysisType::ThreePrime);

                                                if needs_coverage {
                                                    // Helper to process one base position.
                                                    // Only called for positions within [window_start, window_end)
                                                    // so each base is counted by exactly one window.
                                                    let process_pos = |p: u64, state: &mut WorkerState| {
                                                        use crate::analysis::index::Hits;
                                                        match dense.get_hits(p) {
                                                            Hits::None => {}
                                                            Hits::Single(g_idx) => {
                                                                let g_idx = g_idx as usize;
                                                                state.overlaps_found += 1;
                                                                let gene = &genes_ref[g_idx];
                                                                let entry = state.get_or_init(g_idx);
                                                                if let Some(s_5p) = gene.bin_map.get_spliced_5p(p) {
                                                                    entry.add(s_5p as usize, gene.total_len as usize);
                                                                }
                                                            }
                                                            Hits::Multi(indices) => {
                                                                for &g_idx in indices {
                                                                    let g_idx = g_idx as usize;
                                                                    state.overlaps_found += 1;
                                                                    let gene = &genes_ref[g_idx];
                                                                    let entry = state.get_or_init(g_idx);
                                                                    if let Some(s_5p) = gene.bin_map.get_spliced_5p(p) {
                                                                        entry.add(s_5p as usize, gene.total_len as usize);
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    };

                                                if args.ends {
                                                    let is_reverse = flags.is_reverse_complemented();
                                                    let target_pos = if is_reverse {
                                                        let mut curr = pos as u64;
                                                        for op_res in record.cigar().iter() {
                                                            if let Ok(op) = op_res {
                                                                use noodles::sam::alignment::record::cigar::op::Kind;
                                                                if matches!(
                                                                    op.kind(),
                                                                    Kind::Match
                                                                        | Kind::SequenceMatch
                                                                        | Kind::SequenceMismatch
                                                                        | Kind::Deletion
                                                                        | Kind::Skip
                                                                ) {
                                                                    curr += op.len() as u64;
                                                                }
                                                            }
                                                        }
                                                        curr.saturating_sub(1)
                                                    } else {
                                                        pos as u64
                                                    };
                                                    process_pos(target_pos, &mut state);
                                                } else {
                                                    // Full CIGAR walk. Clamp each M-block to
                                                    // [window_start, window_end) so reads that
                                                    // span a window boundary are counted once
                                                    // per base (by whichever window owns that
                                                    // position) rather than once per window.
                                                    let mut curr = pos as u64;
                                                    for op_res in record.cigar().iter() {
                                                        if curr >= win_end { break; }
                                                        if let Ok(op) = op_res {
                                                            use noodles::sam::alignment::record::cigar::op::Kind;
                                                            match op.kind() {
                                                                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                                                                    let block_end = curr + op.len() as u64;
                                                                    // clamp to [win_start, win_end)
                                                                    let p_lo = curr.max(win_start);
                                                                    let p_hi = block_end.min(win_end);
                                                                    for p in p_lo..p_hi {
                                                                        process_pos(p, &mut state);
                                                                    }
                                                                    curr = block_end;
                                                                }
                                                                Kind::Deletion | Kind::Skip => {
                                                                    curr += op.len() as u64;
                                                                }
                                                                _ => {}
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    Err(e) => {
                                            eprintln!(
                                                "    [ERROR] Query record error in {}:{}-{}: {:?}",
                                                chrom, window.start, window.end, e
                                            );
                                            return;
                                        }
                                    }
                                }
                            }
                            Err(e) => {
                                println!(
                                    "    [ERROR] Query failed for region {}:{}-{}: {}",
                                    chrom,
                                    window.start + 1,
                                    window.end,
                                    e
                                );
                            }
                        }
                    });
                    pb.inc(1);
                    state
                })
                .reduce(
                    move || WorkerState::new_with_capacity(n_genes, max_3p_dist, qc_sample_size),
                    WorkerState::merge,
                );

            // Apply merged results to global Gene arrays (via populated list, not full scan).
            for &g_idx in &total_state.populated {
                let g = g_idx as usize;
                if let Some(local) = &total_state.gene_counts[g] {
                    index_proto.genes[g]
                        .apply_local_updates(&local.diff_3p, &local.diff_percentile);
                }
            }
            pb.finish_with_message("Parallel Scan Finished.");

            println!(
                "  - Summary: {} exonic base-level overlaps processed.",
                total_state.overlaps_found
            );
            println!(
                "  - Total records seen by query iterator: {}",
                total_state.records_seen
            );
            println!("    - Failed unmapped:  {}", total_state.fail_unmapped);
            println!("    - Failed secondary: {}", total_state.fail_secondary);
            println!("    - Failed QC:        {}", total_state.fail_qc);
            println!(
                "    - Failed MapQ < {}: {}",
                args.mapq, total_state.fail_mapq
            );
            let passed_filters = total_state.records_seen
                - total_state.fail_unmapped
                - total_state.fail_secondary
                - total_state.fail_qc
                - total_state.fail_mapq;
            println!("    - Passed filters:   {}", passed_filters);
            // Aggregate results for this sample
            if args.analysis.contains(&AnalysisType::GeneBody) {
                let (classic_raw, classic_percent) = aggregate_rseqc_classic(
                    &index_proto,
                    args.ends,
                    args.min_support as u32,
                    args.coverage_weighted,
                );
                classic_all.insert(sample_name.to_string(), classic_raw);
                classic_percent_all.insert(sample_name.to_string(), classic_percent);

                if args.stratify_length && !args.no_plot {
                    println!("  - Length-stratified gene body coverage:");
                    let strat = aggregate_rseqc_stratified(
                        &index_proto,
                        args.min_support as u32,
                        args.coverage_weighted,
                    );
                    if !strat.is_empty() {
                        let strat_path = format!(
                            "{}.{}.geneBodyCoverage.stratified.svg",
                            args.output, sample_name
                        );
                        generate_stratified_gene_body_plot(&strat, &strat_path)?;
                        println!(
                            "  - Length-stratified plot written to: {}",
                            strat_path
                        );
                    }
                }
            }


            if args.analysis.contains(&AnalysisType::ThreePrime) {
                let stats = aggregate_genes(
                    &index_proto,
                    args.normalization_bp,
                    args.min_support,
                    args.max_3p_dist,
                    args.ends,
                );
                dist_3p_all.insert(sample_name.to_string(), stats.dist_3p_means);
                percentile_all.insert(sample_name.to_string(), stats.percentile_means);
            }

            if args.analysis.contains(&AnalysisType::Qc) {
                let mut qc = total_state.qc.summary.clone();
                qc.aligned_qc_reads = total_state.aligned_qc_reads;
                qc.mtdna_reads = total_state.mtdna_reads;
                qc.rdna_reads = total_state.rdna_reads;
                let qc_summary_path = format!("{}.{}.rna_qc.txt", args.output, sample_name);
                qc.write_summary_file(&qc_summary_path)?;
                println!("  - RNA-seq QC summary written to: {}", qc_summary_path);
                println!(
                    "    - mtDNA reads: {} ({:.2}%)",
                    qc.mtdna_reads,
                    qc.mtdna_fraction() * 100.0
                );
                println!(
                    "    - rDNA reads:  {} ({:.2}%)",
                    qc.rdna_reads,
                    qc.rdna_fraction() * 100.0
                );
                println!(
                    "    - Inferred experiment: {} (FR {:.3}, RF {:.3}, other {:.3})",
                    qc.inferred_strandness(),
                    qc.fr_fraction(),
                    qc.rf_fraction(),
                    qc.other_fraction()
                );
                println!(
                    "    - Inner distance mean: {} from {} sampled pairs (target {})",
                    qc.inner_distance_mean()
                        .map(|v| format!("{:.2}", v))
                        .unwrap_or_else(|| "NA".to_string()),
                    qc.informative_pairs,
                    args.qc_sample_size
                );

                let qc_hist_path = format!("{}.{}.inner_distance.tsv", args.output, sample_name);
                qc.write_inner_distance_histogram(&qc_hist_path)?;
                println!("  - Inner distance histogram written to: {}", qc_hist_path);
                inner_distance_all.insert(
                    sample_name.to_string(),
                    qc.clipped_inner_distance_series(-200, 200),
                );
            }

            if args.analysis.contains(&AnalysisType::Distribution) {
                let dist_path = format!("{}.{}.read_distribution.txt", args.output, sample_name);
                let mut f_dist = File::create(&dist_path)?;
                use std::io::Write;

                let total_assigned: u64 = total_state
                    .read_dist_counts
                    .iter()
                    .filter(|(&r, _)| r != RegionType::Intergenic)
                    .map(|(_, &c)| c)
                    .sum();

                writeln!(f_dist, "{:<30}{}", "Total Reads", passed_filters)?;
                writeln!(f_dist, "{:<30}{}", "Total Tags", total_state.total_tags)?;
                writeln!(f_dist, "{:<30}{}", "Total Assigned Tags", total_assigned)?;
                writeln!(
                    f_dist,
                    "====================================================================="
                )?;
                writeln!(
                    f_dist,
                    "{:<20}{:<20}{:<20}{:<20}",
                    "Group", "Total_bases", "Tag_count", "Tags/Kb"
                )?;

                let groups = [
                    RegionType::CdsExon,
                    RegionType::Utr5Exon,
                    RegionType::Utr3Exon,
                    RegionType::Exon,
                    RegionType::Intron,
                    RegionType::TssUp1kb,
                    RegionType::TssUp5kb,
                    RegionType::TssUp10kb,
                    RegionType::TesDown1kb,
                    RegionType::TesDown5kb,
                    RegionType::TesDown10kb,
                ];

                for group in groups {
                    let size = feature_index
                        .as_ref()
                        .and_then(|index| index.feature_sizes.get(&group).cloned())
                        .unwrap_or(0);
                    let count = total_state
                        .read_dist_counts
                        .get(&group)
                        .cloned()
                        .unwrap_or(0);
                    let tags_per_kb = if size > 0 {
                        (count as f64 * 1000.0) / (size as f64)
                    } else {
                        0.0
                    };
                    writeln!(
                        f_dist,
                        "{:<20}{:<20}{:<20}{:<18.2}",
                        group.to_string(),
                        size,
                        count,
                        tags_per_kb
                    )?;
                }
                writeln!(
                    f_dist,
                    "====================================================================="
                )?;
                println!("  - Read distribution report written to: {}", dist_path);
            }

            // RESET index for next sample to avoid cross-contamination
            index_proto.reset_coverage();
        } else {
            println!(
                "  - No BAI Index found. Sequential streaming path not updated for DenseMap yet."
            );
        }
    }

    // Global Plotting & Output
    if args.analysis.contains(&AnalysisType::GeneBody) && !classic_all.is_empty() {
        let classic_path = format!("{}.geneBodyCoverage.txt", args.output);
        write_classic_wide_format(&classic_path, &classic_all)?;
        println!(
            "\nMultiQC-compatible raw aggregate written to: {}",
            classic_path
        );

        let classic_pct_path = format!("{}.geneBodyCoverage.percent.txt", args.output);
        write_classic_wide_format(&classic_pct_path, &classic_percent_all)?;
        println!(
            "Normalized percentage aggregate written to: {}",
            classic_pct_path
        );

        if !args.no_plot {
            generate_gene_body_plot(
                &classic_percent_all,
                &format!("{}.geneBodyCoverage.svg", args.output),
            )?;
            println!(
                "Summary plot written to: {}.geneBodyCoverage.svg",
                args.output
            );
        }
    }

    if args.analysis.contains(&AnalysisType::ThreePrime) && !dist_3p_all.is_empty() {
        if !args.no_plot {
            generate_multi_3p_dist_plot(&dist_3p_all, &format!("{}.3p_dist.svg", args.output))?;
            println!(
                "3' Distance profile written to: {}.3p_dist.svg",
                args.output
            );
        }
    }

    if args.analysis.contains(&AnalysisType::Qc) && !inner_distance_all.is_empty() {
        if !args.no_plot {
            generate_multi_inner_distance_plot(
                &inner_distance_all,
                &format!("{}.inner_distance.svg", args.output),
                -200,
            )?;
            println!(
                "Inner distance profile written to: {}.inner_distance.svg",
                args.output
            );
        }
    }

    println!("\nAll tasks completed in {:.2?}s.", start_time.elapsed());
    Ok(())
}
