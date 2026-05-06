use crate::models::{Exon, Gene, Transcript};
use crate::io::text::open_maybe_gz;
use anyhow::{anyhow, Context, Result};
use regex::Regex;
use std::collections::HashMap;
use std::io::BufRead;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdResolutionStrategy {
    GtfExplicit,
    BedRegex,
    BedDelimiter,
    Fallback(String),
}

impl std::fmt::Display for IdResolutionStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::GtfExplicit => write!(f, "GTF-explicit"),
            Self::BedRegex => write!(f, "BED-regex"),
            Self::BedDelimiter => write!(f, "BED-delimiter"),
            Self::Fallback(reason) => write!(f, "Fallback ({})", reason),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ParsedTranscript {
    pub gene_id: String,
    pub transcript_id: String,
    pub gene_name: Option<String>,
    pub biotype: Option<String>,
    pub chrom: String,
    pub strand: char,
    pub exons: Vec<Exon>,
    pub cds_start: Option<u64>,
    pub cds_end: Option<u64>,
    pub strategy: IdResolutionStrategy,
}

/// Controls which representative isoform is chosen per gene in standard (non-transcript-centric) mode.
///
/// The choice affects the 100-bin gene body coverage plot because the transcript's
/// annotated length determines how bases are binned: selecting a long isoform with an
/// un-expressed 3\' UTR shifts the aggregate profile toward the 5\' end.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum IsoformSelect {
    /// Longest spliced transcript (classic RSeQC behaviour; default).
    #[default]
    Longest,
    /// Shortest qualifying transcript. Avoids long un-expressed 3\' UTRs that
    /// artificially shift the gene-body coverage plot toward the 5\' end.
    Shortest,
    /// Transcript closest to the median length among qualifying isoforms.
    /// A reasonable compromise between completeness and avoiding UTR artefacts.
    Median,
}

pub struct AnnotationConfig {
    pub gene_id_delimiter: Option<char>,
    pub gene_id_regex: Option<String>,
    pub biotype_filter: Option<String>,
    pub three_prime_cluster_window: u64,
    pub min_transcript_length: u64,
    pub max_3p_dist: usize,
    pub three_prime_bin_size: usize,
    pub transcript_centric: bool,
    /// If true, only transcripts on the '+' strand are loaded.
    pub plus_strand_only: bool,
    /// Which isoform to use as representative when a gene has multiple eligible transcripts.
    pub isoform_select: IsoformSelect,
    /// Deprecated compatibility flag. 3' sanity filtering is applied during
    /// coverage aggregation based on observed terminal coverage, not while
    /// loading transcripts, so gene-body coverage can still use the selected
    /// representative isoform.
    pub strict_cluster: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationFormat {
    Auto,
    Bed12,
    Gtf,
}

#[derive(Clone)]
pub struct LoadedAnnotation {
    pub genes: Vec<Gene>,
    pub transcripts: Vec<ParsedTranscript>,
}

pub fn load_genes(
    path: &str,
    format: AnnotationFormat,
    config: &AnnotationConfig,
) -> Result<Vec<Gene>> {
    Ok(load_annotation(path, format, config)?.genes)
}

pub fn load_annotation(
    path: &str,
    format: AnnotationFormat,
    config: &AnnotationConfig,
) -> Result<LoadedAnnotation> {
    let reader = open_maybe_gz(path)
        .with_context(|| format!("Failed to open annotation file: {}", path))?;

    let resolved_format = match format {
        AnnotationFormat::Auto => {
            if path.ends_with(".gtf") || path.ends_with(".gtf.gz") {
                AnnotationFormat::Gtf
            } else {
                AnnotationFormat::Bed12
            }
        }
        f => f,
    };

    let transcripts = match resolved_format {
        AnnotationFormat::Gtf => load_transcripts_from_gtf(reader, config)?,
        AnnotationFormat::Bed12 => load_transcripts_from_bed12(reader, config)?,
        AnnotationFormat::Auto => unreachable!(),
    };

    let genes = build_genes_from_transcripts(&transcripts, config);
    Ok(LoadedAnnotation { genes, transcripts })
}

fn load_transcripts_from_gtf(
    reader: impl BufRead,
    config: &AnnotationConfig,
) -> Result<Vec<ParsedTranscript>> {
    let mut tx_map: HashMap<String, ParsedTranscript> = HashMap::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 {
            continue;
        }

        let feature = fields[2];
        if feature != "exon" && feature != "CDS" {
            continue;
        }

        let chrom = fields[0].to_string();
        let start_1based: u64 = fields[3]
            .parse()
            .map_err(|_| anyhow!("Invalid start at line {}", lineno + 1))?;
        let end_1based: u64 = fields[4]
            .parse()
            .map_err(|_| anyhow!("Invalid end at line {}", lineno + 1))?;
        let strand = fields[6].chars().next().unwrap_or('+');

        // Apply strand filter before doing any further parsing.
        if config.plus_strand_only && strand != '+' {
            continue;
        }

        let attrs = parse_gtf_attributes(fields[8]);
        let gene_id = match attrs.get("gene_id") {
            Some(id) => id.clone(),
            None => {
                if lineno < 1000 {
                    eprintln!("Warning: Missing gene_id at line {}. Skipping.", lineno + 1);
                }
                continue;
            }
        };
        let transcript_id = match attrs.get("transcript_id") {
            Some(id) => id.clone(),
            None => {
                continue;
            }
        };
        let gene_name = attrs.get("gene_name").cloned();

        let biotype = attrs
            .get("transcript_type")
            .or(attrs.get("gene_type"))
            .or(attrs.get("gene_biotype"))
            .or(attrs.get("transcript_biotype"))
            .cloned();

        // Apply biotype filter early if possible (Case-insensitive)
        if let Some(ref target) = config.biotype_filter {
            let targets: Vec<&str> = target.split(',').map(|s| s.trim()).collect();
            if let Some(ref b) = biotype {
                if !targets.iter().any(|&t| t.eq_ignore_ascii_case(b)) {
                    continue;
                }
            } else {
                continue;
            }
        }

        let entry = tx_map
            .entry(transcript_id.clone())
            .or_insert_with(|| ParsedTranscript {
                gene_id,
                transcript_id: transcript_id.clone(),
                gene_name,
                biotype,
                chrom: chrom.clone(),
                strand,
                exons: Vec::new(),
                cds_start: None,
                cds_end: None,
                strategy: IdResolutionStrategy::GtfExplicit,
            });

        let start_0based = start_1based - 1;
        let end_0based = end_1based;

        if feature == "exon" {
            entry.exons.push(Exon {
                start: start_0based,
                end: end_0based,
            });
        } else if feature == "CDS" {
            entry.cds_start = Some(
                entry
                    .cds_start
                    .map_or(start_0based, |s| s.min(start_0based)),
            );
            entry.cds_end = Some(entry.cds_end.map_or(end_0based, |e| e.max(end_0based)));
        }
    }

    let mut transcripts: Vec<ParsedTranscript> = tx_map.into_values().collect();
    for tx in &mut transcripts {
        tx.exons.sort_by_key(|e| e.start);
    }

    Ok(transcripts)
}

fn load_transcripts_from_bed12(
    reader: impl BufRead,
    config: &AnnotationConfig,
) -> Result<Vec<ParsedTranscript>> {
    let regex = config
        .gene_id_regex
        .as_ref()
        .map(|s| Regex::new(s))
        .transpose()?;

    let mut transcripts = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let chrom = fields[0].to_string();
        let start = fields[1].parse::<u64>()?;
        let name = fields[3].to_string();
        let strand = fields[5].chars().next().unwrap_or('+');

        if config.plus_strand_only && strand != '+' {
            continue;
        }
        let thick_start = fields[6].parse::<u64>()?;
        let thick_end = fields[7].parse::<u64>()?;
        let block_count = fields[9].parse::<usize>()?;
        let block_sizes: Vec<u64> = fields[10]
            .trim_end_matches(',')
            .split(',')
            .map(|s| s.parse().unwrap())
            .collect();
        let block_starts: Vec<u64> = fields[11]
            .trim_end_matches(',')
            .split(',')
            .map(|s| s.parse().unwrap())
            .collect();

        // Resolve Gene ID with strategy tracking
        let (gene_id, strategy) = if let Some(re) = &regex {
            if let Some(caps) = re.captures(&name) {
                let id = caps
                    .get(1)
                    .map(|m| m.as_str())
                    .unwrap_or(caps.get(0).unwrap().as_str())
                    .to_string();
                (id, IdResolutionStrategy::BedRegex)
            } else {
                (
                    name.clone(),
                    IdResolutionStrategy::Fallback("Regex no match".to_string()),
                )
            }
        } else if let Some(delim) = config.gene_id_delimiter {
            let id = name.split(delim).next().unwrap_or(&name).to_string();
            (id, IdResolutionStrategy::BedDelimiter)
        } else {
            (
                name.clone(),
                IdResolutionStrategy::Fallback("Transcript-as-gene".to_string()),
            )
        };

        let mut exons = Vec::with_capacity(block_count);
        for i in 0..block_count {
            let e_start = start + block_starts[i];
            let e_end = e_start + block_sizes[i];
            exons.push(Exon {
                start: e_start,
                end: e_end,
            });
        }

        transcripts.push(ParsedTranscript {
            gene_id,
            transcript_id: name,
            gene_name: None,
            biotype: None, // BED12 doesn't typically have biotype
            chrom,
            strand,
            exons,
            cds_start: if thick_start < thick_end {
                Some(thick_start)
            } else {
                None
            },
            cds_end: if thick_start < thick_end {
                Some(thick_end)
            } else {
                None
            },
            strategy,
        });
    }

    Ok(transcripts)
}

fn build_genes_from_transcripts(
    parsed: &[ParsedTranscript],
    config: &AnnotationConfig,
) -> Vec<Gene> {
    let mut by_gene: HashMap<String, Vec<ParsedTranscript>> = HashMap::new();
    for tx in parsed {
        by_gene
            .entry(tx.gene_id.clone())
            .or_default()
            .push(tx.clone());
    }

    let mut genes = Vec::new();
    let mut strategy_counts: HashMap<String, usize> = HashMap::new();
    let mut divergent_count = 0;

    for (_gene_id, tx_meta) in by_gene {
        if tx_meta.is_empty() {
            continue;
        }

        // Track strategy for this gene (use the first one found, they should be consistent)
        let strategy = tx_meta[0].strategy.clone();
        *strategy_counts.entry(strategy.to_string()).or_default() += 1;

        let mut eligible = Vec::new();
        let mut ends = Vec::new();

        for meta in tx_meta {
            let tx = Transcript::new(
                meta.transcript_id.clone(),
                meta.chrom.clone(),
                meta.strand,
                meta.biotype.clone(),
                meta.exons.clone(),
                meta.cds_start,
                meta.cds_end,
            );
            if tx.total_length >= config.min_transcript_length {
                ends.push(tx.three_prime_end());
                eligible.push((meta, tx));
            }
        }

        if eligible.is_empty() {
            continue;
        }

        if config.transcript_centric {
            // TRANSCRIPT-CENTRIC: Add all eligible transcripts as separate Gene objects
            for (meta, tx) in eligible {
                genes.push(Gene::new(
                    meta.gene_id,
                    meta.gene_name,
                    tx,
                    config.max_3p_dist,
                    config.three_prime_bin_size,
                ));
            }
        } else {
            // STANDARD: Track divergent 3' annotations for reporting, then select
            // one representative isoform. Coverage is still collected for the
            // selected transcript; 3' sanity filtering happens later from the
            // observed terminal coverage window.
            ends.sort_unstable();
            let min_3p = ends[0];
            let max_3p = ends[ends.len() - 1];

            if max_3p - min_3p > config.three_prime_cluster_window {
                divergent_count += 1;
            }

            // Sort eligible by spliced length for median/shortest/longest selection.
            eligible.sort_by_key(|(_, tx)| tx.total_length);
            let n = eligible.len();

            let pick = match config.isoform_select {
                IsoformSelect::Longest => n - 1,
                IsoformSelect::Shortest => 0,
                IsoformSelect::Median => n / 2,
            };
            let (target_meta, target_tx) = eligible.swap_remove(pick);

            genes.push(Gene::new(
                target_meta.gene_id,
                target_meta.gene_name,
                target_tx,
                config.max_3p_dist,
                config.three_prime_bin_size,
            ));
        }
    }

    println!("Gene ID Resolution Summary:");
    for (strat, count) in strategy_counts {
        println!("  - {}: {} genes", strat, count);
    }
    if divergent_count > 0 {
        println!(
            "  - Note: Selected representative isoforms for {} genes with divergent 3' ends (>{}bp).",
            divergent_count, config.three_prime_cluster_window
        );
        if config.strict_cluster {
            println!(
                "  - Note: --strict-cluster is handled by observed 3' coverage during aggregation."
            );
        }
    }

    genes
}

fn parse_gtf_attributes(attr: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for field in attr.split(';') {
        let field = field.trim();
        if field.is_empty() {
            continue;
        }

        // Handle both space and = separators (GTF vs GFF-style)
        let (key, val_part) = if let Some(pos) = field.find(' ') {
            (field[..pos].trim(), field[pos..].trim())
        } else if let Some(pos) = field.find('=') {
            (field[..pos].trim(), field[pos..].trim())
        } else {
            continue;
        };

        let val = val_part.trim_matches('"');
        if !key.is_empty() {
            map.insert(key.to_string(), val.to_string());
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::BufReader;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_temp_file(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "{}", content).unwrap();
        file
    }

    #[test]
    fn test_gtf_parsing() {
        let gtf_content = r#"
chr1	custom	exon	101	200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	custom	exon	301	400	.	+	.	gene_id "G1"; transcript_id "T1";
"#;
        let file = create_temp_file(gtf_content);
        let path = file.path().to_str().unwrap();
        let f = File::open(path).unwrap();
        let reader = BufReader::new(f);
        let config = AnnotationConfig {
            gene_id_delimiter: None,
            gene_id_regex: None,
            biotype_filter: None,
            three_prime_cluster_window: 50,
            min_transcript_length: 100,
            max_3p_dist: 15000,
            three_prime_bin_size: 50,
            transcript_centric: false,
            plus_strand_only: false,
            isoform_select: IsoformSelect::Longest,
            strict_cluster: false,
        };
        let txs = load_transcripts_from_gtf(reader, &config).unwrap();

        assert_eq!(txs.len(), 1);
        let tx = &txs[0];
        assert_eq!(tx.gene_id, "G1");
        assert_eq!(tx.transcript_id, "T1");
        assert_eq!(tx.chrom, "chr1");
        assert_eq!(tx.strand, '+');
        assert_eq!(tx.exons.len(), 2);
        // 1-based [101, 200] -> 0-based [100, 200)
        assert_eq!(tx.exons[0].start, 100);
        assert_eq!(tx.exons[0].end, 200);
    }

    #[test]
    fn test_cross_format_consistency() {
        // Same transcript represented in both formats
        // Use G1|T1 so delimiter '|' extracts G1 as gene_id
        let bed_content = "chr1	100	400	G1|T1	0	+	100	400	0	2	100,100	0,200";
        let gtf_content = r#"
chr1	custom	exon	101	200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	custom	exon	301	400	.	+	.	gene_id "G1"; transcript_id "T1";
"#;
        let bed_file = create_temp_file(bed_content);
        let gtf_file = create_temp_file(gtf_content);

        let config = AnnotationConfig {
            gene_id_delimiter: Some('|'),
            gene_id_regex: None,
            biotype_filter: None,
            three_prime_cluster_window: 50,
            min_transcript_length: 100,
            max_3p_dist: 15000,
            three_prime_bin_size: 50,
            transcript_centric: false,
            plus_strand_only: false,
            isoform_select: IsoformSelect::Longest,
            strict_cluster: false,
        };

        let genes_bed = load_genes(
            bed_file.path().to_str().unwrap(),
            AnnotationFormat::Bed12,
            &config,
        )
        .unwrap();
        let genes_gtf = load_genes(
            gtf_file.path().to_str().unwrap(),
            AnnotationFormat::Gtf,
            &config,
        )
        .unwrap();

        assert_eq!(genes_bed.len(), genes_gtf.len());
        assert_eq!(genes_bed[0].id, genes_gtf[0].id);
        assert_eq!(
            genes_bed[0].representative.exons,
            genes_gtf[0].representative.exons
        );
    }

    #[test]
    fn test_gzip_loading() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let gtf_content = r#"
chr1	custom	exon	101	200	.	+	.	gene_id "G1"; transcript_id "T1";
"#;
        let mut file = NamedTempFile::new().unwrap();
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(gtf_content.as_bytes()).unwrap();
        let compressed = encoder.finish().unwrap();
        file.write_all(&compressed).unwrap();

        let config = AnnotationConfig {
            gene_id_delimiter: None,
            gene_id_regex: None,
            biotype_filter: None,
            three_prime_cluster_window: 50,
            min_transcript_length: 100,
            max_3p_dist: 15000,
            three_prime_bin_size: 50,
            transcript_centric: false,
            plus_strand_only: false,
            isoform_select: IsoformSelect::Longest,
            strict_cluster: false,
        };

        // Trick load_genes by using a path that ends in .gz
        // We'll rename the temp file
        let path = file.path().to_str().unwrap();
        let gz_path = format!("{}.gtf.gz", path);
        std::fs::copy(path, &gz_path).unwrap();

        let genes = load_genes(&gz_path, AnnotationFormat::Auto, &config).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].id, "G1");

        std::fs::remove_file(gz_path).unwrap();
    }

    #[test]
    fn test_multi_member_gzip_loading() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        fn gzip_member(content: &str) -> Vec<u8> {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(content.as_bytes()).unwrap();
            encoder.finish().unwrap()
        }

        let member1 = gzip_member(
            r#"chr1	custom	exon	101	200	.	+	.	gene_id "G1"; transcript_id "T1";
"#,
        );
        let member2 = gzip_member(
            r#"chr1	custom	exon	301	400	.	+	.	gene_id "G1"; transcript_id "T1";
"#,
        );

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(&member1).unwrap();
        file.write_all(&member2).unwrap();
        let path = file.path().to_str().unwrap();
        let gz_path = format!("{}.gtf.gz", path);
        std::fs::copy(path, &gz_path).unwrap();

        let config = AnnotationConfig {
            gene_id_delimiter: None,
            gene_id_regex: None,
            biotype_filter: None,
            three_prime_cluster_window: 50,
            min_transcript_length: 100,
            max_3p_dist: 15000,
            three_prime_bin_size: 50,
            transcript_centric: false,
            plus_strand_only: false,
            isoform_select: IsoformSelect::Longest,
            strict_cluster: false,
        };

        let loaded = load_annotation(&gz_path, AnnotationFormat::Auto, &config).unwrap();
        assert_eq!(loaded.transcripts.len(), 1);
        assert_eq!(loaded.transcripts[0].exons.len(), 2);

        std::fs::remove_file(gz_path).unwrap();
    }
}
