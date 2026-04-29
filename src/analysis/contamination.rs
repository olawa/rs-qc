use crate::models::normalize_chrom;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContaminantInterval {
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromContaminantIndex {
    pub intervals: Arc<[ContaminantInterval]>,
}

#[derive(Debug, Default, Clone)]
pub struct ContaminantCursor {
    idx: usize,
}

#[derive(Debug, Default, Clone)]
pub struct ContaminantIndex {
    pub chroms: HashMap<String, Arc<ChromContaminantIndex>>,
}

impl ContaminantIndex {
    pub fn from_bed(path: &str) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("Failed to open rDNA BED file: {}", path))?;
        let reader = BufReader::new(file);
        let mut by_chrom: HashMap<String, Vec<ContaminantInterval>> = HashMap::new();

        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                anyhow::bail!("Invalid BED line {} in {}", line_no + 1, path);
            }

            let chrom = normalize_chrom(fields[0]);
            let start: u64 = fields[1]
                .parse()
                .with_context(|| format!("Invalid BED start at line {}", line_no + 1))?;
            let end: u64 = fields[2]
                .parse()
                .with_context(|| format!("Invalid BED end at line {}", line_no + 1))?;

            if start < end {
                by_chrom
                    .entry(chrom)
                    .or_default()
                    .push(ContaminantInterval { start, end });
            }
        }

        let chroms = by_chrom
            .into_iter()
            .filter_map(|(chrom, intervals)| {
                let merged = merge_intervals(intervals);
                if merged.is_empty() {
                    None
                } else {
                    Some((
                        chrom,
                        Arc::new(ChromContaminantIndex {
                            intervals: merged.into(),
                        }),
                    ))
                }
            })
            .collect();

        Ok(Self { chroms })
    }
}

impl ChromContaminantIndex {
    pub fn cursor_at(&self, pos: u64) -> ContaminantCursor {
        ContaminantCursor {
            idx: self.intervals.partition_point(|iv| iv.end <= pos),
        }
    }

    pub fn contains(&self, pos: u64, cursor: &mut ContaminantCursor) -> bool {
        while cursor.idx < self.intervals.len() && self.intervals[cursor.idx].end <= pos {
            cursor.idx += 1;
        }

        self.intervals
            .get(cursor.idx)
            .is_some_and(|iv| iv.start <= pos && pos < iv.end)
    }
}

fn merge_intervals(mut intervals: Vec<ContaminantInterval>) -> Vec<ContaminantInterval> {
    intervals.sort_by_key(|iv| (iv.start, iv.end));
    let mut merged: Vec<ContaminantInterval> = Vec::new();

    for interval in intervals {
        if let Some(last) = merged.last_mut() {
            if interval.start <= last.end {
                last.end = last.end.max(interval.end);
                continue;
            }
        }
        merged.push(interval);
    }

    merged
}
