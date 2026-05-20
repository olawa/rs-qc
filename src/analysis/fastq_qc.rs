use crate::analysis::report::write_summary_json;
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use kuva::prelude::*;
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;
use std::process::{Child, ChildStdout, Command, Stdio};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;

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
    pub threads: usize,
    pub sample_size: usize,
    pub kmer_size: usize,
    pub top_n: usize,
    pub phred_offset: u8,
    pub no_kmers: bool,
    pub paired: bool,
    pub length_bin_size: usize,
    pub use_pigz: bool,
    pub pigz_threads: usize,
    pub batch_size: usize,
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

#[derive(Clone, Debug, Default)]
pub struct BasePositionMetrics {
    pub count: u64,
    pub qual_sum: u64,
    pub bases: [u64; 5],
    pub qual_hist: BTreeMap<u8, u64>,
    pub adapter_hits: u64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct BasePositionSummary {
    pub count: u64,
    pub mean_quality: f64,
    pub bases: [u64; 5],
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
    pub per_base: Vec<BasePositionSummary>,
}

#[derive(Debug, Serialize)]
pub struct FastqLengthBin {
    pub bin_start: usize,
    pub bin_end: usize,
    pub reads: u64,
    pub bases: u64,
    pub cumulative_bases_ge_start: u64,
}

#[derive(Clone, Debug)]
struct FastqRecordOwned {
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    sampled: bool,
}

#[derive(Debug)]
struct FastqBatch {
    records: Vec<FastqRecordOwned>,
}

#[derive(Debug)]
struct PairedFastqBatch {
    left: Vec<FastqRecordOwned>,
    right: Vec<FastqRecordOwned>,
}

#[derive(Debug)]
enum FastqWorkItem {
    Single(FastqBatch),
    Paired(PairedFastqBatch),
}

impl FastqQcMetrics {
    #[cfg(test)]
    fn observe(&mut self, seq: &[u8], qual: &[u8], config: &FastqQcConfig) -> Result<()> {
        self.observe_record(
            &FastqRecordOwned {
                name: Vec::new(),
                seq: seq.to_vec(),
                qual: qual.to_vec(),
                sampled: true,
            },
            config,
        )
    }

    fn observe_record(&mut self, record: &FastqRecordOwned, config: &FastqQcConfig) -> Result<()> {
        let seq = record.seq.as_slice();
        let qual = record.qual.as_slice();
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

        if record.sampled {
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

    fn merge(mut self, other: Self) -> Self {
        self.total_reads += other.total_reads;
        self.total_bases += other.total_bases;

        if self.min_len == 0 || (other.min_len > 0 && other.min_len < self.min_len) {
            self.min_len = other.min_len;
        }
        self.max_len = self.max_len.max(other.max_len);

        for (k, v) in other.length_hist {
            *self.length_hist.entry(k).or_insert(0) += v;
        }
        for (k, v) in other.gc_hist {
            *self.gc_hist.entry(k).or_insert(0) += v;
        }
        for (k, v) in other.mean_quality_hist {
            *self.mean_quality_hist.entry(k).or_insert(0) += v;
        }
        for (k, v) in other.adapter_hits {
            *self.adapter_hits.entry(k).or_insert(0) += v;
        }
        for (k, v) in other.overrepresented_bytes {
            *self.overrepresented_bytes.entry(k).or_insert(0) += v;
        }
        for (k, v) in other.kmers {
            *self.kmers.entry(k).or_insert(0) += v;
        }

        self.duplicate_sample_reads += other.duplicate_sample_reads;
        self.duplicate_sample_unique += other.duplicate_sample_unique;
        self.poly_a_reads += other.poly_a_reads;
        self.poly_g_reads += other.poly_g_reads;
        self.n_reads += other.n_reads;
        self.paired_reads_checked += other.paired_reads_checked;
        self.paired_name_mismatches += other.paired_name_mismatches;

        if self.per_base.len() < other.per_base.len() {
            self.per_base
                .resize_with(other.per_base.len(), BasePositionMetrics::default);
        }
        for (idx, other_pos) in other.per_base.into_iter().enumerate() {
            let pos = &mut self.per_base[idx];
            pos.count += other_pos.count;
            pos.qual_sum += other_pos.qual_sum;
            for i in 0..pos.bases.len() {
                pos.bases[i] += other_pos.bases[i];
            }
            for (q, c) in other_pos.qual_hist {
                *pos.qual_hist.entry(q).or_insert(0) += c;
            }
            pos.adapter_hits += other_pos.adapter_hits;
        }

        self
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
            per_base: self
                .per_base
                .iter()
                .map(|p| {
                    let count = p.count.max(1);
                    BasePositionSummary {
                        count: p.count,
                        mean_quality: p.qual_sum as f64 / count as f64,
                        bases: p.bases,
                        adapter_hits: p.adapter_hits,
                    }
                })
                .collect(),
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
        metrics = metrics.merge(scan_paired_fastq(
            &config.inputs[0],
            &config.inputs[1],
            config,
        )?);
    } else {
        let mut observed_reads = 0_u64;
        for input in &config.inputs {
            let sample_metrics = scan_fastq(input, config, &mut observed_reads)
                .with_context(|| format!("failed while scanning FASTQ input {input}"))?;
            metrics = metrics.merge(sample_metrics);
        }
    }

    write_outputs(config, &metrics)?;
    Ok(metrics)
}

fn scan_paired_fastq(
    left_path: &str,
    right_path: &str,
    config: &FastqQcConfig,
) -> Result<FastqQcMetrics> {
    scan_fastq_pairs(left_path, right_path, config)
}

fn scan_fastq(
    path: &str,
    config: &FastqQcConfig,
    observed_reads: &mut u64,
) -> Result<FastqQcMetrics> {
    let (metrics, next_observed_reads) = scan_fastq_single(path, config, *observed_reads)?;
    *observed_reads = next_observed_reads;
    Ok(metrics)
}

fn scan_fastq_single(
    path: &str,
    config: &FastqQcConfig,
    observed_reads: u64,
) -> Result<(FastqQcMetrics, u64)> {
    let reader = BufReader::with_capacity(
        4 << 20, // 4MB buffer for maximum I/O performance
        open_fastq_reader(path, config.use_pigz, config.pigz_threads)?,
    );
    let producer = FastqBatchProducer::new_single(reader, config, observed_reads)?;
    producer.run()
}

fn scan_fastq_pairs(
    left_path: &str,
    right_path: &str,
    config: &FastqQcConfig,
) -> Result<FastqQcMetrics> {
    let left = BufReader::with_capacity(
        4 << 20, // 4MB buffer
        open_fastq_reader(left_path, config.use_pigz, config.pigz_threads)?,
    );
    let right = BufReader::with_capacity(
        4 << 20, // 4MB buffer
        open_fastq_reader(right_path, config.use_pigz, config.pigz_threads)?,
    );
    let producer = FastqBatchProducer::new_paired(left, right, config, 0)?;
    let (metrics, _) = producer.run()?;
    Ok(metrics)
}

struct FastqBatchProducer {
    left: BufReader<Box<dyn Read>>,
    right: Option<BufReader<Box<dyn Read>>>,
    config: FastqQcConfig,
    observed_reads: u64,
    tx: Option<mpsc::SyncSender<FastqWorkItem>>,
    rx: Arc<Mutex<mpsc::Receiver<FastqWorkItem>>>,
    start_time: std::time::Instant,
    last_progress: u64,
}

impl FastqBatchProducer {
    fn new_single(
        reader: BufReader<Box<dyn Read>>,
        config: &FastqQcConfig,
        observed_reads: u64,
    ) -> Result<Self> {
        Self::new(reader, None, config, observed_reads)
    }

    fn new_paired(
        left: BufReader<Box<dyn Read>>,
        right: BufReader<Box<dyn Read>>,
        config: &FastqQcConfig,
        observed_reads: u64,
    ) -> Result<Self> {
        Self::new(left, Some(right), config, observed_reads)
    }

    fn new(
        left: BufReader<Box<dyn Read>>,
        right: Option<BufReader<Box<dyn Read>>>,
        config: &FastqQcConfig,
        observed_reads: u64,
    ) -> Result<Self> {
        let worker_threads = worker_thread_count(config);
        let queue_capacity = worker_threads.saturating_mul(2).max(1);
        let (tx, rx) = mpsc::sync_channel(queue_capacity);
        let start_time = std::time::Instant::now();
        Ok(Self {
            left,
            right,
            config: config.clone(),
            observed_reads,
            tx: Some(tx),
            rx: Arc::new(Mutex::new(rx)),
            start_time,
            last_progress: 0,
        })
    }

    fn run(mut self) -> Result<(FastqQcMetrics, u64)> {
        let worker_count = worker_thread_count(&self.config);
        let worker_config = Arc::new(self.config.clone());
        let rx = Arc::clone(&self.rx);
        let mut handles = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let rx = Arc::clone(&rx);
            let worker_config = Arc::clone(&worker_config);
            handles.push(thread::spawn(move || worker_loop(rx, worker_config)));
        }

        if self.right.is_some() {
            self.produce_paired_batches()?;
        } else {
            self.produce_single_batches()?;
        }

        drop(self.tx.take());

        let mut merged = FastqQcMetrics::default();
        for handle in handles {
            let worker_metrics = handle
                .join()
                .map_err(|_| anyhow::anyhow!("FASTQ worker thread panicked"))??;
            merged = merged.merge(worker_metrics);
        }
        Ok((merged, self.observed_reads))
    }

    fn produce_single_batches(&mut self) -> Result<()> {
        let mut buffer = FastqRecordBuffer::default();
        let mut batch = FastqBatch {
            records: Vec::with_capacity(self.config.batch_size.max(1)),
        };
        loop {
            match read_record(
                &mut self.left,
                &mut buffer,
                self.config.sample_size,
                &mut self.observed_reads,
                "FASTQ",
            )? {
                Some(record) => {
                    batch.records.push(record);
                    if batch.records.len() >= self.config.batch_size.max(1) {
                        self.tx
                            .as_ref()
                            .unwrap()
                            .send(FastqWorkItem::Single(batch))
                            .map_err(|_| {
                                anyhow::anyhow!("FASTQ workers stopped receiving batches")
                            })?;
                        batch = FastqBatch {
                            records: Vec::with_capacity(self.config.batch_size.max(1)),
                        };
                    }

                    // Periodic stderr progress update
                    let report_threshold = 1_000_000;
                    if self.observed_reads >= self.last_progress + report_threshold {
                        self.last_progress = (self.observed_reads / report_threshold) * report_threshold;
                        let elapsed = self.start_time.elapsed().as_secs_f64();
                        let speed = if elapsed > 0.0 { self.observed_reads as f64 / elapsed } else { 0.0 };
                        eprintln!(
                            "  - Processed {} reads ({:.0} reads/s)...",
                            self.observed_reads, speed
                        );
                    }
                }
                None => break,
            }
        }
        if !batch.records.is_empty() {
            self.tx
                .as_ref()
                .unwrap()
                .send(FastqWorkItem::Single(batch))
                .map_err(|_| anyhow::anyhow!("FASTQ workers stopped receiving batches"))?;
        }
        Ok(())
    }

    fn produce_paired_batches(&mut self) -> Result<()> {
        let right = self.right.as_mut().unwrap();
        let mut left_buffer = FastqRecordBuffer::default();
        let mut right_buffer = FastqRecordBuffer::default();
        let mut batch = PairedFastqBatch {
            left: Vec::with_capacity(self.config.batch_size.max(1)),
            right: Vec::with_capacity(self.config.batch_size.max(1)),
        };
        loop {
            let left_record = read_record(
                &mut self.left,
                &mut left_buffer,
                self.config.sample_size,
                &mut self.observed_reads,
                "left FASTQ",
            )?;
            let right_record = read_record(
                right,
                &mut right_buffer,
                self.config.sample_size,
                &mut self.observed_reads,
                "right FASTQ",
            )?;
            match (left_record, right_record) {
                (None, None) => break,
                (Some(_), None) | (None, Some(_)) => {
                    bail!("paired FASTQ files have different record counts")
                }
                (Some(left), Some(right)) => {
                    batch.left.push(left);
                    batch.right.push(right);
                    if batch.left.len() >= self.config.batch_size.max(1) {
                        self.tx
                            .as_ref()
                            .unwrap()
                            .send(FastqWorkItem::Paired(batch))
                            .map_err(|_| {
                                anyhow::anyhow!("FASTQ workers stopped receiving batches")
                            })?;
                        batch = PairedFastqBatch {
                            left: Vec::with_capacity(self.config.batch_size.max(1)),
                            right: Vec::with_capacity(self.config.batch_size.max(1)),
                        };
                    }

                    // Periodic stderr progress update
                    let report_threshold = 1_000_000;
                    if self.observed_reads >= self.last_progress + report_threshold {
                        self.last_progress = (self.observed_reads / report_threshold) * report_threshold;
                        let elapsed = self.start_time.elapsed().as_secs_f64();
                        let speed = if elapsed > 0.0 { self.observed_reads as f64 / elapsed } else { 0.0 };
                        eprintln!(
                            "  - Processed {} reads ({:.0} reads/s)...",
                            self.observed_reads, speed
                        );
                    }
                }
            }
        }

        if !batch.left.is_empty() {
            self.tx
                .as_ref()
                .unwrap()
                .send(FastqWorkItem::Paired(batch))
                .map_err(|_| anyhow::anyhow!("FASTQ workers stopped receiving batches"))?;
        }
        Ok(())
    }
}

struct FastqRecordBuffer {
    name: Vec<u8>,
    seq: Vec<u8>,
    plus: Vec<u8>,
    qual: Vec<u8>,
}

impl Default for FastqRecordBuffer {
    fn default() -> Self {
        Self {
            name: Vec::new(),
            seq: Vec::new(),
            plus: Vec::new(),
            qual: Vec::new(),
        }
    }
}

fn worker_loop(
    rx: Arc<Mutex<mpsc::Receiver<FastqWorkItem>>>,
    config: Arc<FastqQcConfig>,
) -> Result<FastqQcMetrics> {
    let mut metrics = FastqQcMetrics::default();
    loop {
        let work = {
            let guard = rx.lock().unwrap();
            guard.recv()
        };
        match work {
            Ok(FastqWorkItem::Single(batch)) => {
                process_single_batch(&mut metrics, batch, &config)?;
            }
            Ok(FastqWorkItem::Paired(batch)) => {
                process_paired_batch(&mut metrics, batch, &config)?;
            }
            Err(_) => break,
        }
    }
    Ok(metrics)
}

fn process_single_batch(
    metrics: &mut FastqQcMetrics,
    batch: FastqBatch,
    config: &FastqQcConfig,
) -> Result<()> {
    for record in batch.records {
        metrics.observe_record(&record, config)?;
    }
    Ok(())
}

fn process_paired_batch(
    metrics: &mut FastqQcMetrics,
    batch: PairedFastqBatch,
    config: &FastqQcConfig,
) -> Result<()> {
    for (left, right) in batch.left.into_iter().zip(batch.right.into_iter()) {
        metrics.paired_reads_checked += 1;
        if normalized_read_name(&left.name) != normalized_read_name(&right.name) {
            metrics.paired_name_mismatches += 1;
        }
        metrics.observe_record(&left, config)?;
        metrics.observe_record(&right, config)?;
    }
    Ok(())
}

fn read_record<R: BufRead + ?Sized>(
    reader: &mut R,
    record: &mut FastqRecordBuffer,
    sample_size: usize,
    observed_reads: &mut u64,
    path: &str,
) -> Result<Option<FastqRecordOwned>> {
    record.name.clear();
    if read_fastq_line(reader, &mut record.name)? == 0 {
        return Ok(None);
    }
    record.seq.clear();
    record.plus.clear();
    record.qual.clear();
    if read_fastq_line(reader, &mut record.seq)? == 0
        || read_fastq_line(reader, &mut record.plus)? == 0
        || read_fastq_line(reader, &mut record.qual)? == 0
    {
        bail!("truncated FASTQ record in {path}");
    }

    if !record.name.starts_with(b"@") {
        bail!("invalid FASTQ record in {path}: header does not start with @");
    }
    if !record.plus.starts_with(b"+") {
        bail!("invalid FASTQ record in {path}: separator does not start with +");
    }

    let sampled = *observed_reads < sample_size as u64;
    *observed_reads += 1;

    Ok(Some(FastqRecordOwned {
        name: std::mem::take(&mut record.name),
        seq: std::mem::take(&mut record.seq),
        qual: std::mem::take(&mut record.qual),
        sampled,
    }))
}

fn read_fastq_line<R: BufRead + ?Sized>(reader: &mut R, buf: &mut Vec<u8>) -> Result<usize> {
    buf.clear();
    let n = reader.read_until(b'\n', buf)?;
    if n == 0 {
        return Ok(0);
    }
    while matches!(buf.last(), Some(b'\n' | b'\r')) {
        buf.pop();
    }
    Ok(n)
}

fn open_fastq_reader(path: &str, use_pigz: bool, pigz_threads: usize) -> Result<Box<dyn Read>> {
    if use_pigz && path.ends_with(".gz") {
        match open_pigz_reader(path, pigz_threads) {
            Ok(reader) => return Ok(Box::new(reader)),
            Err(err) => {
                eprintln!(
                    "pigz unavailable for {path}: {err}; falling back to internal gzip reader"
                );
            }
        }
    }
    let file = File::open(path).with_context(|| format!("could not open {path}"))?;
    if path.ends_with(".gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

fn open_pigz_reader(path: &str, pigz_threads: usize) -> Result<PigzReader> {
    let mut child = Command::new("pigz")
        .arg("-dc")
        .arg("-p")
        .arg(pigz_threads.max(1).to_string())
        .arg(path)
        .stdout(Stdio::piped())
        .stdin(Stdio::null())
        .spawn()?;
    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| anyhow::anyhow!("pigz did not provide stdout"))?;
    Ok(PigzReader { child, stdout })
}

struct PigzReader {
    child: Child,
    stdout: ChildStdout,
}

impl Read for PigzReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.stdout.read(buf)
    }
}

impl Drop for PigzReader {
    fn drop(&mut self) {
        let _ = self.child.wait();
    }
}

fn worker_thread_count(config: &FastqQcConfig) -> usize {
    let threads = config.threads.max(1);
    if config.use_pigz {
        threads.saturating_sub(config.pigz_threads.max(1)).max(1)
    } else {
        threads
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
    Ok(())
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

fn normalized_read_name(name: &[u8]) -> String {
    let name = String::from_utf8_lossy(name);
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
            threads: 1,
            sample_size: 10,
            kmer_size: 3,
            top_n: 10,
            phred_offset: 33,
            no_kmers: false,
            paired: false,
            length_bin_size: 1000,
            use_pigz: false,
            pigz_threads: 1,
            batch_size: 16,
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
