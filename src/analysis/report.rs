use anyhow::{bail, Context, Result};
use serde::Serialize;
use serde_json::Value;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Serialize)]
pub struct ReportDocument {
    pub title: String,
    pub generated_by: String,
    pub sections: Vec<ReportSection>,
}

#[derive(Debug, Serialize)]
pub struct ReportSection {
    pub module: String,
    pub sample: String,
    pub source: String,
    pub metrics: Value,
}

#[derive(Debug, Serialize)]
pub struct SummaryEnvelope<T> {
    pub module: String,
    pub sample: String,
    pub metrics: T,
}

pub fn write_summary_json<T: Serialize>(
    path: &str,
    module: &str,
    sample: &str,
    metrics: &T,
) -> Result<()> {
    let envelope = SummaryEnvelope {
        module: module.to_string(),
        sample: sample.to_string(),
        metrics,
    };
    let file = fs::File::create(path).with_context(|| format!("could not create {path}"))?;
    serde_json::to_writer_pretty(file, &envelope)?;
    Ok(())
}

pub fn build_document(source_paths: &[PathBuf]) -> Result<ReportDocument> {
    if source_paths.is_empty() {
        bail!("no summary JSON files were found for the report");
    }

    let mut sections = Vec::new();
    for path in source_paths {
        let raw = fs::read_to_string(path)
            .with_context(|| format!("could not read summary JSON {}", path.display()))?;
        let value: Value = serde_json::from_str(&raw)
            .with_context(|| format!("could not parse summary JSON {}", path.display()))?;

        let (module, sample, metrics) = if let Some(obj) = value.as_object() {
            let module = obj
                .get("module")
                .and_then(|v| v.as_str())
                .unwrap_or_else(|| infer_module_from_path(path))
                .to_string();
            let sample = obj
                .get("sample")
                .and_then(|v| v.as_str())
                .unwrap_or_else(|| infer_sample_from_path(path))
                .to_string();
            let metrics = obj.get("metrics").cloned().unwrap_or(value.clone());
            (module, sample, metrics)
        } else {
            (
                infer_module_from_path(path).to_string(),
                infer_sample_from_path(path).to_string(),
                value.clone(),
            )
        };

        sections.push(ReportSection {
            module,
            sample,
            source: path.display().to_string(),
            metrics,
        });
    }

    Ok(ReportDocument {
        title: "rs-qc summary report".to_string(),
        generated_by: env!("CARGO_PKG_VERSION").to_string(),
        sections,
    })
}

pub fn write_report(document: &ReportDocument, output_prefix: &str) -> Result<()> {
    let json_path = format!("{}.summary.json", output_prefix);
    let html_path = format!("{}.report.html", output_prefix);

    fs::write(&json_path, serde_json::to_string_pretty(document)?)
        .with_context(|| format!("could not write {json_path}"))?;
    fs::write(&html_path, render_html(document))
        .with_context(|| format!("could not write {html_path}"))?;
    Ok(())
}

pub fn resolve_input_files(inputs: &[String]) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for input in inputs {
        let path = Path::new(input);
        if path.exists() && path.extension().and_then(|e| e.to_str()) == Some("json") {
            files.push(path.to_path_buf());
            continue;
        }

        let candidates = [
            format!("{input}.fastq.summary.json"),
            format!("{input}.align.summary.json"),
            format!("{input}.rna.summary.json"),
            format!("{input}.summary.json"),
        ];

        let mut found = false;
        for candidate in candidates {
            let candidate_path = PathBuf::from(&candidate);
            if candidate_path.exists() {
                files.push(candidate_path);
                found = true;
            }
        }

        if !found {
            let parent = path.parent().unwrap_or_else(|| Path::new("."));
            let prefix = path.file_name().and_then(|s| s.to_str()).unwrap_or(input);
            for entry in fs::read_dir(parent)
                .with_context(|| format!("could not read directory {}", parent.display()))?
            {
                let entry = entry?;
                let candidate_path = entry.path();
                let Some(name) = candidate_path.file_name().and_then(|s| s.to_str()) else {
                    continue;
                };
                if name.starts_with(prefix) && name.ends_with(".summary.json") {
                    files.push(candidate_path);
                    found = true;
                }
            }
        }

        if !found {
            bail!(
                "could not find a summary JSON for input `{input}`; pass a summary .json file or an output prefix"
            );
        }
    }

    Ok(files)
}

fn render_html(document: &ReportDocument) -> String {
    let module_counts = overview_module_counts(document);
    let mut out = String::new();
    out.push_str("<!doctype html><html><head><meta charset=\"utf-8\">");
    out.push_str("<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">");
    out.push_str("<title>rs-qc report</title>");
    out.push_str(
        "<style>\
        :root{color-scheme:light}\
        body{font-family:system-ui,-apple-system,BlinkMacSystemFont,\"Segoe UI\",sans-serif; margin:0; padding:32px; background:linear-gradient(180deg,#f8faff 0%,#eef2f8 100%); color:#142033;}\
        .wrap{max-width:1200px; margin:0 auto;}\
        h1{margin:0 0 8px 0; font-size:2.15rem;}\
        .sub{color:#5a6475; margin:0 0 24px 0;}\
        .overview{display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:14px; margin:0 0 24px 0;}\
        .stat{background:#fff; border:1px solid #e2e7f0; border-radius:14px; padding:14px 16px; box-shadow:0 1px 3px rgba(20,32,51,.06);}\
        .stat .label{display:block; color:#64748b; font-size:.82rem; text-transform:uppercase; letter-spacing:.04em; margin-bottom:6px;}\
        .stat .value{font-size:1.4rem; font-weight:700;}\
        .section{background:#fff; border:1px solid #e2e7f0; border-radius:14px; padding:18px 20px; margin:0 0 18px 0; box-shadow:0 1px 3px rgba(20,32,51,.06);}\
        .section-head{display:flex; justify-content:space-between; gap:12px; align-items:flex-start; flex-wrap:wrap;}\
        .badge{display:inline-flex; align-items:center; gap:6px; border-radius:999px; background:#eef4ff; color:#274b9f; padding:4px 10px; font-size:.78rem; font-weight:700; text-transform:uppercase; letter-spacing:.05em;}\
        .meta{display:flex; gap:12px; flex-wrap:wrap; color:#5a6475; font-size:.92rem; margin:10px 0 12px 0;}\
        .links{display:flex; gap:10px; flex-wrap:wrap; margin:10px 0 14px 0;}\
        .links a{display:inline-block; text-decoration:none; color:#1d4ed8; background:#eff6ff; border:1px solid #dbeafe; border-radius:10px; padding:6px 10px; font-size:.9rem;}\
        table{width:100%; border-collapse:collapse; font-size:.94rem;}\
        td{padding:6px 8px; border-top:1px solid #edf1f7; vertical-align:top;}\
        td.key{font-weight:600; width:34%; color:#243047;}\
        details{margin-top:14px;}\
        pre{white-space:pre-wrap; word-break:break-word; background:#0f172a; color:#e5eefc; border-radius:12px; padding:14px; overflow:auto;}\
        </style>",
    );
    out.push_str("</head><body><div class=\"wrap\">");
    out.push_str(&format!(
        "<h1>{}</h1><p class=\"sub\">Generated by rs-qc {}</p>",
        escape_html(&document.title),
        escape_html(&document.generated_by)
    ));
    out.push_str("<div class=\"overview\">");
    out.push_str(&format!(
        "<div class=\"stat\"><span class=\"label\">Modules</span><span class=\"value\">{}</span></div>",
        document.sections.len()
    ));
    for (module, count) in module_counts {
        out.push_str(&format!(
            "<div class=\"stat\"><span class=\"label\">{}</span><span class=\"value\">{}</span></div>",
            escape_html(&module),
            count
        ));
    }
    out.push_str("</div>");

    for section in &document.sections {
        out.push_str("<section class=\"section\">");
        out.push_str("<div class=\"section-head\">");
        out.push_str(&format!(
            "<div><h2 style=\"margin:0\">{}</h2><div class=\"meta\"><span class=\"badge\">{}</span><span>{}</span></div></div>",
            escape_html(&section.sample),
            escape_html(&section.module),
            escape_html(&section.source)
        ));
        out.push_str("</div>");
        let links = artifact_links(section);
        if !links.is_empty() {
            out.push_str("<div class=\"links\">");
            for (label, href) in links {
                out.push_str(&format!(
                    "<a href=\"{}\">{}</a>",
                    escape_html(&href),
                    escape_html(&label)
                ));
            }
            out.push_str("</div>");
        }
        out.push_str("<table><tbody>");
        for (key, value) in flatten_summary(&section.metrics, "") {
            out.push_str(&format!(
                "<tr><td class=\"key\">{}</td><td>{}</td></tr>",
                escape_html(&key),
                escape_html(&value)
            ));
        }
        out.push_str("</tbody></table>");
        out.push_str("<details><summary>Raw JSON</summary><pre>");
        out.push_str(&escape_html(
            &serde_json::to_string_pretty(&section.metrics).unwrap_or_else(|_| "{}".to_string()),
        ));
        out.push_str("</pre></details></section>");
    }

    out.push_str("</div></body></html>");
    out
}

fn overview_module_counts(document: &ReportDocument) -> Vec<(String, usize)> {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for section in &document.sections {
        *counts.entry(section.module.clone()).or_insert(0) += 1;
    }
    counts.into_iter().collect()
}

fn artifact_links(section: &ReportSection) -> Vec<(String, String)> {
    let source = Path::new(&section.source);
    let Some(stem) = source.file_stem().and_then(|s| s.to_str()) else {
        return Vec::new();
    };
    let dir = source.parent().unwrap_or_else(|| Path::new("."));
    let candidates = match section.module.as_str() {
        "fastq" => vec![
            format!("{stem}.fastq.summary.json"),
            format!("{stem}.fastq.summary.txt"),
            format!("{stem}.fastq.per_base.tsv"),
            format!("{stem}.fastq.length_distribution.tsv"),
            format!("{stem}.fastq.gc_distribution.tsv"),
            format!("{stem}.fastq.mean_quality_distribution.tsv"),
            format!("{stem}.fastq.overrepresented.tsv"),
            format!("{stem}.fastq.kmers.tsv"),
        ],
        "align" => vec![
            format!("{stem}.align.summary.json"),
            format!("{stem}.align.summary.txt"),
            format!("{stem}.align.mapq.tsv"),
            format!("{stem}.align.read_length.tsv"),
            format!("{stem}.align.insert_size.tsv"),
            format!("{stem}.align.cigar.tsv"),
            format!("{stem}.align.contigs.tsv"),
            format!("{stem}.align.de_accuracy.tsv"),
        ],
        "rna" => vec![
            format!("{stem}.rna.summary.json"),
            format!("{stem}.rna_qc.txt"),
            format!("{stem}.inner_distance.tsv"),
            format!("{stem}.geneBodyCoverage.txt"),
            format!("{stem}.geneBodyCoverage.svg"),
        ],
        _ => Vec::new(),
    };

    candidates
        .into_iter()
        .filter_map(|name| {
            let path = dir.join(&name);
            if path.exists() {
                Some((name, path.display().to_string()))
            } else {
                None
            }
        })
        .collect()
}

fn flatten_summary(value: &Value, prefix: &str) -> Vec<(String, String)> {
    let mut rows = Vec::new();
    match value {
        Value::Object(map) => {
            for (key, child) in map {
                let next_prefix = if prefix.is_empty() {
                    key.to_string()
                } else {
                    format!("{prefix}.{key}")
                };
                rows.extend(flatten_summary(child, &next_prefix));
            }
        }
        Value::Array(items) => {
            let rendered = if items.iter().all(is_scalar) {
                items
                    .iter()
                    .map(render_scalar)
                    .collect::<Vec<_>>()
                    .join(", ")
            } else {
                format!("{} items", items.len())
            };
            rows.push((prefix.to_string(), rendered));
        }
        _ => rows.push((prefix.to_string(), render_scalar(value))),
    }
    rows
}

fn is_scalar(value: &Value) -> bool {
    matches!(
        value,
        Value::Null | Value::Bool(_) | Value::Number(_) | Value::String(_)
    )
}

fn render_scalar(value: &Value) -> String {
    match value {
        Value::Null => "NA".to_string(),
        Value::Bool(v) => v.to_string(),
        Value::Number(v) => v.to_string(),
        Value::String(v) => v.clone(),
        other => serde_json::to_string(other).unwrap_or_else(|_| "<unrenderable>".to_string()),
    }
}

fn escape_html(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn infer_module_from_path(path: &Path) -> &str {
    let name = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
    if name.contains(".fastq.summary.json") {
        "fastq"
    } else if name.contains(".align.summary.json") {
        "align"
    } else if name.contains(".rna.summary.json") {
        "rna"
    } else if name.contains(".rna") {
        "rna"
    } else {
        "unknown"
    }
}

fn infer_sample_from_path(path: &Path) -> &str {
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("sample")
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::env;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn resolve_input_files_finds_prefixes_and_report_builds_document() {
        let cwd = env::current_dir().expect("current dir");
        let tmp = tempdir().expect("tempdir");
        env::set_current_dir(tmp.path()).expect("set temp cwd");

        let fastq_path = tmp.path().join("sample.fastq.summary.json");
        let align_path = tmp.path().join("sample.align.summary.json");
        let rna_path = tmp.path().join("sample.rna.summary.json");

        write_summary_json(
            fastq_path.to_str().unwrap(),
            "fastq",
            "sample",
            &json!({"reads": 10, "mean_read_length": 100.0}),
        )
        .expect("write fastq summary");
        write_summary_json(
            align_path.to_str().unwrap(),
            "align",
            "sample",
            &json!({"total_records": 25, "mapped_records": 20}),
        )
        .expect("write align summary");
        write_summary_json(
            rna_path.to_str().unwrap(),
            "rna",
            "sample",
            &json!({"aligned_qc_reads": 5, "mtdna_reads": 1}),
        )
        .expect("write rna summary");

        let files = resolve_input_files(&["sample".to_string()]).expect("resolve prefix");
        assert_eq!(files.len(), 3);

        let document = build_document(&files).expect("build document");
        assert_eq!(document.sections.len(), 3);
        assert!(document.sections.iter().any(|s| s.module == "fastq"));
        assert!(document.sections.iter().any(|s| s.module == "align"));
        assert!(document.sections.iter().any(|s| s.module == "rna"));

        write_report(&document, "report_out").expect("write report");
        assert!(tmp.path().join("report_out.summary.json").exists());
        assert!(tmp.path().join("report_out.report.html").exists());

        env::set_current_dir(cwd).expect("restore cwd");
        fs::remove_dir_all(tmp.path()).ok();
    }
}
