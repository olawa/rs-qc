use crate::models::Exon;
use anyhow::Result;
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::HashMap;
use std::fs::File;
use std::num::NonZeroUsize;

pub enum QcMolecule {
    Single {
        chrom: String,
        blocks: Vec<Exon>,
        is_reverse: bool,
    },
    Paired {
        chrom: String,
        read1_blocks: Vec<Exon>,
        read2_blocks: Vec<Exon>,
        read1_is_reverse: bool,
        read2_is_reverse: bool,
    },
}

pub struct BamFragmentIterator {
    reader: bam::io::Reader<noodles::bgzf::MultithreadedReader<File>>,
    header: sam::Header,
    pending_mates: HashMap<String, bam::Record>,
    mapq_threshold: u8,
    orphans_count: u64,
    discordant_count: u64,
    records_processed: u64,
}

impl BamFragmentIterator {
    pub fn new(path: &str, mapq_threshold: u8, worker_count: usize) -> Result<Self> {
        let file = File::open(path)?;
        let worker_count = NonZeroUsize::new(worker_count).unwrap_or(NonZeroUsize::new(1).unwrap());
        let mt_reader = noodles::bgzf::MultithreadedReader::with_worker_count(worker_count, file);
        let mut reader = bam::io::Reader::from(mt_reader);
        let header = reader.read_header()?;

        Ok(Self {
            reader,
            header,
            pending_mates: HashMap::new(),
            mapq_threshold,
            orphans_count: 0,
            discordant_count: 0,
            records_processed: 0,
        })
    }

    pub fn records_processed(&self) -> u64 {
        self.records_processed
    }

    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    pub fn orphans_count(&self) -> u64 {
        self.orphans_count + self.pending_mates.len() as u64
    }

    pub fn discordant_count(&self) -> u64 {
        self.discordant_count
    }
}

impl Iterator for BamFragmentIterator {
    type Item = Result<QcMolecule>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = bam::Record::default();

        while self.reader.read_record(&mut record).is_ok() {
            self.records_processed += 1;
            let flags = record.flags();
            if flags.is_secondary()
                || flags.is_supplementary()
                || flags.is_unmapped()
                || flags.is_qc_fail()
                || flags.is_duplicate()
            {
                continue;
            }

            let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(0);
            if mapq < self.mapq_threshold {
                continue;
            }

            let qname = record
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).to_string())
                .unwrap_or_default();

            if !flags.is_segmented() {
                let blocks = get_aligned_blocks(&record);
                let chrom = get_chrom(&self.header, &record);
                return Some(Ok(QcMolecule::Single {
                    chrom,
                    blocks,
                    is_reverse: flags.is_reverse_complemented(),
                }));
            }

            if let Some(mate) = self.pending_mates.remove(&qname) {
                let id_rec = record
                    .reference_sequence_id()
                    .and_then(|res| res.ok())
                    .map(|id| usize::from(id));
                let id_mate = mate
                    .reference_sequence_id()
                    .and_then(|res| res.ok())
                    .map(|id| usize::from(id));

                if id_rec == id_mate && id_rec.is_some() && is_sane_orientation(&record, &mate) {
                    let chrom = get_chrom(&self.header, &record);
                    let record_is_first = flags.is_first_segment();
                    let mate_flags = mate.flags();

                    let (read1_blocks, read1_is_reverse, read2_blocks, read2_is_reverse) =
                        if record_is_first
                            || (!mate_flags.is_first_segment() && !mate_flags.is_last_segment())
                        {
                            (
                                get_aligned_blocks(&record),
                                flags.is_reverse_complemented(),
                                get_aligned_blocks(&mate),
                                mate_flags.is_reverse_complemented(),
                            )
                        } else {
                            (
                                get_aligned_blocks(&mate),
                                mate_flags.is_reverse_complemented(),
                                get_aligned_blocks(&record),
                                flags.is_reverse_complemented(),
                            )
                        };

                    return Some(Ok(QcMolecule::Paired {
                        chrom,
                        read1_blocks,
                        read2_blocks,
                        read1_is_reverse,
                        read2_is_reverse,
                    }));
                } else {
                    self.discordant_count += 1;
                    continue;
                }
            } else {
                if self.pending_mates.len() < 1_000_000 {
                    self.pending_mates.insert(qname, record.clone());
                } else {
                    static WARNED: std::sync::atomic::AtomicBool =
                        std::sync::atomic::AtomicBool::new(false);
                    if !WARNED.swap(true, std::sync::atomic::Ordering::SeqCst) {
                        eprintln!("WARNING: Fragment buffer exceeded 1M records. Are you using a coordinate-sorted BAM with many orphans? Dropping new orphans to save memory.");
                    }
                    self.orphans_count += 1;
                    let blocks = get_aligned_blocks(&record);
                    let chrom = get_chrom(&self.header, &record);
                    return Some(Ok(QcMolecule::Single {
                        chrom,
                        blocks,
                        is_reverse: flags.is_reverse_complemented(),
                    }));
                }
            }
        }
        None
    }
}

fn is_sane_orientation(r1: &bam::Record, r2: &bam::Record) -> bool {
    let f1 = r1.flags();
    let f2 = r2.flags();
    f1.is_reverse_complemented() != f2.is_reverse_complemented()
}

fn get_aligned_blocks(record: &bam::Record) -> Vec<Exon> {
    let mut blocks = Vec::new();
    let start = match record.alignment_start() {
        Some(Ok(s)) => s.get() as u64 - 1,
        _ => return blocks,
    };

    let mut curr_pos = start;
    let cigar = record.cigar();
    for op_res in cigar.iter() {
        if let Ok(op) = op_res {
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let end = curr_pos + op.len() as u64;
                    blocks.push(Exon {
                        start: curr_pos,
                        end,
                    });
                    curr_pos = end;
                }
                Kind::Deletion | Kind::Skip => {
                    curr_pos += op.len() as u64;
                }
                _ => {}
            }
        }
    }
    blocks
}

fn get_chrom(header: &sam::Header, record: &bam::Record) -> String {
    if let Some(Ok(id)) = record.reference_sequence_id() {
        if let Some((name, _)) = header.reference_sequences().get_index(id) {
            return String::from_utf8_lossy(name.as_ref()).to_string();
        }
    }
    "unknown".to_string()
}
