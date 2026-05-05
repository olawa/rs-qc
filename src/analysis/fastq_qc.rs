use crate::analysis::report::write_summary_json;
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use kuva::prelude::*;
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

const DEFAULT_ADAPTERS: &[&str] = &[
    "AGATCGGAAGAGC",
    "CTGTCTCTTATACACATCT",
    "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    "AATGATACGGCGACCACCGAGATCTACAC",
];

#[derive(Clone, Debug)]
pub struct FastqQcConfig {
    pub inputs: Vec<String>,
    pub output_prefix: String,
    pub sample_size: usize,
    pub kmer_size: usize,
    pub top_n: usize,
    pub phred_offset: u8,
    pub no_kmers: bool,
    pub paired: bool,
    pub length_bin_size: usize,
}

#[derive(Debug, Default)]
pub struct FastqQcMetrics {
    pub total_reads: u64,
    pub total_bases: u64,
    pub min_len: usize,
    pub max_len: usize,
    pub per_base: Vec<BasePositionMetrics>,
    pub length_hist: BTreeMap<usize, u64>,
    pub gc_hist: BTreeMap<u8, u64>,
    pub mean_quality_hist: BTreeMap<u8, u64>,
    pub adapter_hits: HashMap<String, u64>,
    pub overrepresented_bytes: HashMap<Vec<u8>, u64>,
    pub kmers: HashMap<Vec<u8>, u64>,
    pub duplicate_sample_reads: u64,
    pub duplicate_sample_unique: u64,
    pub poly_a_reads: u64,
    pub poly_g_reads: u64,
    pub n_reads: u64,
    pub paired_reads_checked: u64,
    pub paired_name_mismatches: u64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct BasePositionMetrics {
    pub count: u64,
    pub qual_sum: u64,
    pub bases: [u64; 5],
    pub qual_hist: BTreeMap<u8, u64>,
    pub adapter_hits: u64,
}

#[derive(Debug, Serialize)]
pub struct FastqKmerCount {
    pub kmer: String,
    pub count: u64,
}

#[derive(Debug, Serialize)]
pub struct FastqOverrepresented {
    pub sequence: String,
    pub count: u64,
}

#[derive(Debug, Serialize)]
pub struct FastqQcSummary {
    pub total_reads: u64,
    pub total_bases: u64,
    pub min_len: usize,
    pub max_len: usize,
    pub mean_read_length: f64,
    pub read_nx: BTreeMap<u8, usize>,
    pub length_hist: BTreeMap<usize, u64>,
    pub gc_hist: BTreeMap<u8, u64>,
    pub mean_quality_hist: BTreeMap<u8, u64>,
    pub adapter_hits: BTreeMap<String, u64>,
    pub overrepresented: Vec<FastqOverrepresented>,
    pub kmers: Vec<FastqKmerCount>,
    pub duplicate_sample_reads: u64,
    pub duplicate_sample_unique: u64,
    pub duplication_estimate: f64,
    pub poly_a_reads: u64,
    pub poly_g_reads: u64,
    pub n_reads: u64,
    pub paired_reads_checked: u64,
    pub paired_name_mismatches: u64,
    pub per_base: Vec<BasePositionMetrics>,
}

#[derive(Debug, Serialize)]
pub struct FastqLengthBin {
    pub bin_start: usize,
    pub bin_end: usize,
    pub reads: u64,
    pub bases: u64,
    pub cumulative_bases_ge_start: u64,
}

impl FastqQcMetrics {
    fn observe(&mut self, seq: &[u8], qual: &[u8], config: &FastqQcConfig) -> Result<()> {
        if seq.len() != qual.len() {
            bail!(
                "FASTQ sequence and quality lengths differ: {} sequence bases vs {} quality scores",
                seq.len(),
                qual.len()
            );
        }

        self.total_reads += 1;
        self.total_bases += seq.len() as u64;
        if self.total_reads == 1 {
            self.min_len = seq.len();
            self.max_len = seq.len();
        } else {
            self.min_len = self.min_len.min(seq.len());
            self.max_len = self.max_len.max(seq.len());
        }
        *self.length_hist.entry(seq.len()).or_insert(0) += 1;

        let mut gc = 0_u64;
        let mut n_seen = false;
        let mut qual_sum = 0_u64;
        if self.per_base.len() < seq.len() {
            self.per_base
                .resize_with(seq.len(), BasePositionMetrics::default);
        }

        for (i, (&base, &q)) in seq.iter().zip(qual.iter()).enumerate() {
            let phred = q.saturating_sub(config.phred_offset) as u64;
            qual_sum += phred;

            let idx = match base.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => {
                    gc += 1;
                    1
                }
                b'G' => {
                    gc += 1;
                    2
                }
                b'T' | b'U' => 3,
                _ => {
                    n_seen = true;
                    4
                }
            };

            let pos = &mut self.per_base[i];
            pos.count += 1;
            pos.qual_sum += phred;
            pos.bases[idx] += 1;
            *pos.qual_hist.entry(phred as u8).or_insert(0) += 1;
        }

        if n_seen {
            self.n_reads += 1;
        }

        if has_poly_tail(seq, b'A', 12) {
            self.poly_a_reads += 1;
        }
        if has_poly_tail(seq, b'G', 12) {
            self.poly_g_reads += 1;
        }

        let gc_pct = if seq.is_empty() {
            0
        } else {
            ((gc * 100) / seq.len() as u64) as u8
        };
        *self.gc_hist.entry(gc_pct).or_insert(0) += 1;

        let mean_q = if seq.is_empty() {
            0
        } else {
            (qual_sum / seq.len() as u64) as u8
        };
        *self.mean_quality_hist.entry(mean_q).or_insert(0) += 1;

        for adapter in DEFAULT_ADAPTERS {
            if let Some(found_pos) = find_substring(seq, adapter.as_bytes()) {
                *self.adapter_hits.entry((*adapter).to_string()).or_insert(0) += 1;
                if found_pos < self.per_base.len() {
                    self.per_base[found_pos].adapter_hits += 1;
                }
            }
        }

        if self.duplicate_sample_reads < config.sample_size as u64 {
            self.duplicate_sample_reads += 1;

            // Use entry pattern without to_string() where possible
            match self.overrepresented_bytes.get_mut(seq) {
                Some(count) => *count += 1,
                None => {
                    self.duplicate_sample_unique += 1;
                    self.overrepresented_bytes.insert(seq.to_vec(), 1);
                }
            }

            if !config.no_kmers && config.kmer_size > 0 && seq.len() >= config.kmer_size {
                let mut kmer_buf = vec![0u8; config.kmer_size];
                for kmer in seq.windows(config.kmer_size) {
                    let mut valid = true;
                    for (j, &b) in kmer.iter().enumerate() {
                        let upper = b.to_ascii_uppercase();
                        if !matches!(upper, b'A' | b'C' | b'G' | b'T') {
                            valid = false;
                            break;
                        }
                        kmer_buf[j] = upper;
                    }
                    if valid {
                        *self.kmers.entry(kmer_buf.clone()).or_insert(0) += 1;
                    }
                }
            }
        }

        Ok(())
    }

    pub fn mean_read_length(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            self.total_bases as f64 / self.total_reads as f64
        }
    }

    pub fn duplication_estimate(&self) -> f64 {
        if self.duplicate_sample_reads == 0 {
            0.0
        } else {
            1.0 - (self.duplicate_sample_unique as f64 / self.duplicate_sample_reads as f64)
        }
    }

    pub fn read_nx(&self) -> BTreeMap<u8, usize> {
        let mut out = BTreeMap::new();
        if self.total_bases == 0 {
            return out;
        }

        let mut next_n = 5_u8;
        let mut cumulative_bases = 0_u64;
        for (&len, &count) in self.length_hist.iter().rev() {
            cumulative_bases += len as u64 * count;
            while next_n <= 100
                && cumulative_bases.saturating_mul(100)
                    >= self.total_bases.saturating_mul(next_n as u64)
            {
                out.insert(next_n, len);
                next_n += 5;
            }
        }
        out
    }

    pub fn length_bins(&self, bin_size: usize) -> Vec<FastqLengthBin> {
        let bin_size = bin_size.max(1);
        let mut bins: BTreeMap<usize, (u64, u64)> = BTreeMap::new();
        for (&len, &count) in &self.length_hist {
            let bin_start = (len / bin_size) * bin_size;
            let entry = bins.entry(bin_start).or_default();
            entry.0 += count;
            entry.1 += count * len as u64;
        }

        let mut cumulative = 0_u64;
        let mut out = Vec::with_capacity(bins.len());
        for (&bin_start, &(reads, bases)) in bins.iter().rev() {
            cumulative += bases;
            out.push(FastqLengthBin {
                bin_start,
                bin_end: bin_start + bin_size,
                reads,
                bases,
                cumulative_bases_ge_start: cumulative,
            });
        }
        out.reverse();
        out
    }

    pub fn summary(&self) -> FastqQcSummary {
        let mut overrepresented: Vec<_> = self
            .overrepresented_bytes
            .iter()
            .map(|(sequence, &count)| FastqOverrepresented {
                sequence: String::from_utf8_lossy(sequence).to_string(),
                count,
            })
            .collect();
        overrepresented.sort_by(|a, b| {
            b.count
                .cmp(&a.count)
                .then_with(|| a.sequence.cmp(&b.sequence))
        });

        let mut kmers: Vec<_> = self
            .kmers
            .iter()
            .map(|(kmer, &count)| FastqKmerCount {
                kmer: String::from_utf8_lossy(kmer).to_string(),
                count,
            })
            .collect();
        kmers.sort_by(|a, b| b.count.cmp(&a.count).then_with(|| a.kmer.cmp(&b.kmer)));

        FastqQcSummary {
            total_reads: self.total_reads,
            total_bases: self.total_bases,
            min_len: self.min_len,
            max_len: self.max_len,
            mean_read_length: self.mean_read_length(),
            read_nx: self.read_nx(),
            length_hist: self.length_hist.clone(),
            gc_hist: self.gc_hist.clone(),
            mean_quality_hist: self.mean_quality_hist.clone(),
            adapter_hits: self
                .adapter_hits
                .iter()
                .map(|(k, &v)| (k.clone(), v))
                .collect(),
            overrepresented,
            kmers,
            duplicate_sample_reads: self.duplicate_sample_reads,
            duplicate_sample_unique: self.duplicate_sample_unique,
            duplication_estimate: self.duplication_estimate(),
            poly_a_reads: self.poly_a_reads,
            poly_g_reads: self.poly_g_reads,
            n_reads: self.n_reads,
            paired_reads_checked: self.paired_reads_checked,
            paired_name_mismatches: self.paired_name_mismatches,
            per_base: self.per_base.clone(),
        }
    }
}

pub fn run_fastq_qc(config: &FastqQcConfig) -> Result<FastqQcMetrics> {
    if config.inputs.is_empty() {
        bail!("at least one FASTQ input is required");
    }

    let mut metrics = FastqQcMetrics::default();
    if config.paired {
        if config.inputs.len() != 2 {
            bail!("--paired requires exactly two FASTQ inputs");
        }
        scan_paired_fastq(&config.inputs[0], &config.inputs[1], config, &mut metrics)?;
    } else {
        for input in &config.inputs {
            scan_fastq(input, config, &mut metrics)
                .with_context(|| format!("failed while scanning FASTQ input {input}"))?;
        }
    }

    write_outputs(config, &metrics)?;
    Ok(metrics)
}

fn scan_paired_fastq(
    left_path: &str,
    right_path: &str,
    config: &FastqQcConfig,
    metrics: &mut FastqQcMetrics,
) -> Result<()> {
    let mut left = open_fastq(left_path)?;
    let mut right = open_fastq(right_path)?;
    let mut l = FastqRecordBuffer::default();
    let mut r = FastqRecordBuffer::default();

    loop {
        let left_record = read_record(&mut left, &mut l, left_path)?;
        let right_record = read_record(&mut right, &mut r, right_path)?;
        match (left_record, right_record) {
            (false, false) => break,
            (true, false) | (false, true) => {
                bail!("paired FASTQ files have different record counts")
            }
            (true, true) => {
                metrics.paired_reads_checked += 1;
                if normalized_read_name(&l.name) != normalized_read_name(&r.name) {
                    metrics.paired_name_mismatches += 1;
                }
                metrics.observe(
                    l.seq.trim_end().as_bytes(),
                    l.qual.trim_end().as_bytes(),
                    config,
                )?;
                metrics.observe(
                    r.seq.trim_end().as_bytes(),
                    r.qual.trim_end().as_bytes(),
                    config,
                )?;
            }
        }
    }

    Ok(())
}

fn scan_fastq(path: &str, config: &FastqQcConfig, metrics: &mut FastqQcMetrics) -> Result<()> {
    let mut reader = open_fastq(path)?;
    let mut record = FastqRecordBuffer::default();

    while read_record(&mut reader, &mut record, path)? {
        metrics.observe(
            record.seq.trim_end().as_bytes(),
            record.qual.trim_end().as_bytes(),
            config,
        )?;
    }

    Ok(())
}

#[derive(Default)]
struct FastqRecordBuffer {
    name: String,
    seq: String,
    plus: String,
    qual: String,
}

fn read_record<R: BufRead + ?Sized>(
    reader: &mut R,
    record: &mut FastqRecordBuffer,
    path: &str,
) -> Result<bool> {
    record.name.clear();
    if reader.read_line(&mut record.name)? == 0 {
        return Ok(false);
    }
    record.seq.clear();
    record.plus.clear();
    record.qual.clear();
    if reader.read_line(&mut record.seq)? == 0
        || reader.read_line(&mut record.plus)? == 0
        || reader.read_line(&mut record.qual)? == 0
    {
        bail!("truncated FASTQ record in {path}");
    }

    if !record.name.starts_with('@') {
        bail!("invalid FASTQ record in {path}: header does not start with @");
    }
    if !record.plus.starts_with('+') {
        bail!("invalid FASTQ record in {path}: separator does not start with +");
    }

    Ok(true)
}

fn open_fastq(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("could not open {path}"))?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn write_outputs(config: &FastqQcConfig, metrics: &FastqQcMetrics) -> Result<()> {
    let summary = metrics.summary();
    write_summary(config, metrics)?;
    write_per_base(config, metrics)?;
    write_histogram(
        &format!("{}.fastq.length_distribution.tsv", config.output_prefix),
        "length",
        "reads",
        metrics.length_hist.iter().map(|(&k, &v)| (k as u64, v)),
    )?;
    write_length_bins(config, metrics)?;
    write_length_plot(config, metrics)?;
    write_histogram(
        &format!("{}.fastq.gc_distribution.tsv", config.output_prefix),
        "gc_percent",
        "reads",
        metrics.gc_hist.iter().map(|(&k, &v)| (k as u64, v)),
    )?;
    write_histogram(
        &format!(
            "{}.fastq.mean_quality_distribution.tsv",
            config.output_prefix
        ),
        "mean_quality",
        "reads",
        metrics
            .mean_quality_hist
            .iter()
            .map(|(&k, &v)| (k as u64, v)),
    )?;
    write_ranked_sequences(
        &format!("{}.fastq.overrepresented.tsv", config.output_prefix),
        "sequence",
        metrics
            .overrepresented_bytes
            .iter()
            .map(|(sequence, &count)| (std::str::from_utf8(sequence).unwrap_or("N"), count)),
        config.top_n,
    )?;
    write_ranked_sequences(
        &format!("{}.fastq.kmers.tsv", config.output_prefix),
        "kmer",
        metrics
            .kmers
            .iter()
            .map(|(k, &v)| (std::str::from_utf8(k).unwrap_or("N"), v)),
        config.top_n,
    )?;
    write_summary_json(
        &format!("{}.fastq.summary.json", config.output_prefix),
        "fastq",
        &config.output_prefix,
        &summary,
    )?;
    write_per_base_metrics(config, metrics)?;
    Ok(())
}

fn write_per_base_metrics(config: &FastqQcConfig, metrics: &FastqQcMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.fastq.per_base.tsv", config.output_prefix))?;
    writeln!(
        out,
        "pos\tcount\tA\tC\tG\tT\tN\tmean_qual\tmedian_qual\tadapter_hits"
    )?;
    for (i, p) in metrics.per_base.iter().enumerate() {
        let mean_q = if p.count > 0 {
            p.qual_sum as f64 / p.count as f64
        } else {
            0.0
        };
        let median_q = calculate_median(&p.qual_hist);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
            i + 1,
            p.count,
            p.bases[0],
            p.bases[1],
            p.bases[2],
            p.bases[3],
            p.bases[4],
            mean_q,
            median_q,
            p.adapter_hits
        )?;
    }
    Ok(())
}

fn calculate_median(hist: &BTreeMap<u8, u64>) -> u8 {
    let total: u64 = hist.values().sum();
    if total == 0 {
        return 0;
    }
    let mut seen = 0;
    let target = total / 2;
    for (&q, &c) in hist {
        seen += c;
        if seen >= target {
            return q;
        }
    }
    0
}

fn write_summary(config: &FastqQcConfig, metrics: &FastqQcMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.fastq.summary.txt", config.output_prefix))?;
    writeln!(out, "total_reads\t{}", metrics.total_reads)?;
    writeln!(out, "total_bases\t{}", metrics.total_bases)?;
    writeln!(out, "min_read_length\t{}", metrics.min_len)?;
    writeln!(out, "mean_read_length\t{:.2}", metrics.mean_read_length())?;
    writeln!(out, "max_read_length\t{}", metrics.max_len)?;
    for (n, len) in metrics.read_nx() {
        writeln!(out, "read_n{n}\t{len}")?;
    }
    writeln!(
        out,
        "duplication_estimate\t{:.6}",
        metrics.duplication_estimate()
    )?;
    writeln!(
        out,
        "sampled_reads_for_duplicates\t{}",
        metrics.duplicate_sample_reads
    )?;
    writeln!(out, "poly_a_reads\t{}", metrics.poly_a_reads)?;
    writeln!(out, "poly_g_reads\t{}", metrics.poly_g_reads)?;
    writeln!(out, "n_containing_reads\t{}", metrics.n_reads)?;
    writeln!(
        out,
        "paired_reads_checked\t{}",
        metrics.paired_reads_checked
    )?;
    writeln!(
        out,
        "paired_name_mismatches\t{}",
        metrics.paired_name_mismatches
    )?;
    for (adapter, count) in &metrics.adapter_hits {
        writeln!(out, "adapter_hit:{adapter}\t{count}")?;
    }
    Ok(())
}

fn write_length_bins(config: &FastqQcConfig, metrics: &FastqQcMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.fastq.length_bins.tsv", config.output_prefix))?;
    writeln!(
        out,
        "bin_start\tbin_end\treads\tbases\tcumulative_bases_ge_start"
    )?;
    for bin in metrics.length_bins(config.length_bin_size) {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            bin.bin_start, bin.bin_end, bin.reads, bin.bases, bin.cumulative_bases_ge_start
        )?;
    }
    Ok(())
}

fn write_length_plot(config: &FastqQcConfig, metrics: &FastqQcMetrics) -> Result<()> {
    let bins = metrics.length_bins(config.length_bin_size);
    if bins.is_empty() {
        return Ok(());
    }

    let data: Vec<(f64, f64)> = bins
        .iter()
        .map(|bin| (bin.bin_start as f64, bin.reads as f64))
        .collect();
    let line = LinePlot::new()
        .with_data(data)
        .with_legend("reads".to_string())
        .with_line_style(LineStyle::Solid);
    let plots = vec![line.into()];
    let layout = Layout::auto_from_plots(&plots)
        .with_title("FASTQ Read Length Distribution")
        .with_x_label("Read length bin start (bp)")
        .with_y_label("Reads");
    let svg = render_to_svg(plots, layout);
    std::fs::write(
        format!("{}.fastq.length_distribution.svg", config.output_prefix),
        svg,
    )?;
    Ok(())
}

fn write_per_base(config: &FastqQcConfig, metrics: &FastqQcMetrics) -> Result<()> {
    let mut out = File::create(format!("{}.fastq.per_base.tsv", config.output_prefix))?;
    writeln!(
        out,
        "position\tcount\tmean_quality\tA\tC\tG\tT\tN\tA_pct\tC_pct\tG_pct\tT_pct\tN_pct"
    )?;
    for (i, pos) in metrics.per_base.iter().enumerate() {
        let count = pos.count.max(1);
        let mean_q = pos.qual_sum as f64 / count as f64;
        writeln!(
            out,
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            i + 1,
            pos.count,
            mean_q,
            pos.bases[0],
            pos.bases[1],
            pos.bases[2],
            pos.bases[3],
            pos.bases[4],
            pct(pos.bases[0], count),
            pct(pos.bases[1], count),
            pct(pos.bases[2], count),
            pct(pos.bases[3], count),
            pct(pos.bases[4], count)
        )?;
    }
    Ok(())
}

fn write_histogram<I>(path: &str, key: &str, value: &str, rows: I) -> Result<()>
where
    I: IntoIterator<Item = (u64, u64)>,
{
    let mut out = File::create(path)?;
    writeln!(out, "{key}\t{value}")?;
    for (k, v) in rows {
        writeln!(out, "{k}\t{v}")?;
    }
    Ok(())
}

fn write_ranked_sequences<'a, I>(path: &str, label: &str, rows: I, top_n: usize) -> Result<()>
where
    I: IntoIterator<Item = (&'a str, u64)>,
{
    let mut rows: Vec<_> = rows.into_iter().collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));

    let mut out = File::create(path)?;
    writeln!(out, "{label}\tcount")?;
    for (seq, count) in rows.into_iter().take(top_n) {
        writeln!(out, "{seq}\t{count}")?;
    }
    Ok(())
}

fn pct(n: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        (n as f64 * 100.0) / total as f64
    }
}

fn find_substring(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    haystack
        .windows(needle.len())
        .position(|window| window.eq_ignore_ascii_case(needle))
}

fn contains_ascii(haystack: &[u8], needle: &[u8]) -> bool {
    find_substring(haystack, needle).is_some()
}

fn has_poly_tail(seq: &[u8], base: u8, min_run: usize) -> bool {
    let base = base.to_ascii_uppercase();
    let mut run = 0;
    for b in seq.iter().rev() {
        if b.to_ascii_uppercase() == base {
            run += 1;
            if run >= min_run {
                return true;
            }
        } else {
            break;
        }
    }
    false
}

fn normalized_read_name(name: &str) -> String {
    let name = name.trim_start_matches('@').trim();
    let name = name.split_whitespace().next().unwrap_or(name);
    name.trim_end_matches("/1")
        .trim_end_matches("/2")
        .to_string()
}

pub fn sample_name_from_path(path: &str) -> String {
    Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(path)
        .trim_end_matches(".gz")
        .trim_end_matches(".fastq")
        .trim_end_matches(".fq")
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn observes_basic_fastq_metrics() {
        let config = FastqQcConfig {
            inputs: Vec::new(),
            output_prefix: "unused".to_string(),
            sample_size: 10,
            kmer_size: 3,
            top_n: 10,
            phred_offset: 33,
            no_kmers: false,
            paired: false,
            length_bin_size: 1000,
        };
        let mut metrics = FastqQcMetrics::default();
        metrics.observe(b"ACGTNN", b"IIIIII", &config).unwrap();
        metrics.observe(b"ACGTAA", b"!!!!!!", &config).unwrap();

        assert_eq!(metrics.total_reads, 2);
        assert_eq!(metrics.total_bases, 12);
        assert_eq!(metrics.length_hist.get(&6), Some(&2));
        assert_eq!(metrics.n_reads, 1);
        assert_eq!(metrics.per_base[0].bases[0], 2);
        assert_eq!(metrics.per_base[1].bases[1], 2);
        assert_eq!(metrics.mean_quality_hist.get(&40), Some(&1));
        assert_eq!(metrics.mean_quality_hist.get(&0), Some(&1));
    }

    #[test]
    fn computes_long_read_nx_and_length_bins() {
        let mut metrics = FastqQcMetrics::default();
        for len in [100_usize, 200, 700, 1000] {
            metrics.total_reads += 1;
            metrics.total_bases += len as u64;
            *metrics.length_hist.entry(len).or_insert(0) += 1;
        }

        assert_eq!(metrics.read_nx().get(&50), Some(&1000));
        assert_eq!(metrics.read_nx().get(&75), Some(&700));
        let bins = metrics.length_bins(500);
        assert_eq!(bins.len(), 3);
        assert_eq!(bins[0].bin_start, 0);
        assert_eq!(bins[0].reads, 2);
        assert_eq!(bins[2].bin_start, 1000);
        assert_eq!(bins[2].cumulative_bases_ge_start, 1000);
    }
}
