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
            push_unique(&mut files, path.to_path_buf());
            continue;
        }

        let mut found_any = false;
        let candidates = [
            format!("{input}.fastq.summary.json"),
            format!("{input}.align.summary.json"),
            format!("{input}.dna.summary.json"),
            format!("{input}.rna.summary.json"),
            format!("{input}.summary.json"),
        ];

        for candidate in candidates {
            let candidate_path = PathBuf::from(&candidate);
            if candidate_path.exists() {
                push_unique(&mut files, candidate_path);
                found_any = true;
            }
        }

        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
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
                push_unique(&mut files, candidate_path);
                found_any = true;
            }
        }

        if !found_any {
            bail!(
                "could not find a summary JSON for input `{input}`; pass a summary .json file or an output prefix"
            );
        }
    }

    Ok(files)
}

fn push_unique(files: &mut Vec<PathBuf>, candidate: PathBuf) {
    let normalized = fs::canonicalize(&candidate).unwrap_or(candidate);
    if !files
        .iter()
        .any(|existing| fs::canonicalize(existing).ok().as_ref() == Some(&normalized))
    {
        files.push(normalized);
    }
}

fn prune_heavy_fields(val: &mut serde_json::Value) {
    match val {
        serde_json::Value::Object(map) => {
            map.remove("kmers");
            map.remove("overrepresented_bytes");
            map.remove("overrepresented");
            map.remove("per_base");
            for v in map.values_mut() {
                prune_heavy_fields(v);
            }
        }
        serde_json::Value::Array(arr) => {
            for v in arr.iter_mut() {
                prune_heavy_fields(v);
            }
        }
        _ => {}
    }
}

fn render_html(document: &ReportDocument) -> String {
    let mut val = serde_json::to_value(document).unwrap_or(serde_json::Value::Null);
    prune_heavy_fields(&mut val);
    let json_data = serde_json::to_string(&val).unwrap_or_else(|_| "{}".to_string());
    let safe_json_data = json_data.replace("</script>", "<\\/script>");

    let mut out = String::new();
    out.push_str("<!doctype html><html><head><meta charset=\"utf-8\">");
    out.push_str("<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">");
    out.push_str("<title>rs-qc Interactive Dashboard</title>");
    
    // Inject OutFit Google Font and Chart.js from CDN (deferred to prevent offline hanging)
    out.push_str("<link rel=\"preconnect\" href=\"https://fonts.googleapis.com\">");
    out.push_str("<link rel=\"preconnect\" href=\"https://fonts.gstatic.com\" crossorigin>");
    out.push_str("<link href=\"https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=Plus+Jakarta+Sans:wght@300;400;500;600;700&display=swap\" rel=\"stylesheet\">");
    out.push_str("<script src=\"https://cdn.jsdelivr.net/npm/chart.js\" defer></script>");

    // Sleek modern styling
    out.push_str(
        "<style>\
        :root {\
            --primary: #2563eb;\
            --primary-light: #eff6ff;\
            --primary-border: #dbeafe;\
            --success: #10b981;\
            --success-light: #ecfdf5;\
            --warning: #f59e0b;\
            --danger: #ef4444;\
            --bg-dark: #0f172a;\
            --bg-sidebar: #1e293b;\
            --bg-body: #f8fafc;\
            --bg-card: #ffffff;\
            --text-main: #1e293b;\
            --text-muted: #64748b;\
            --border-color: #e2e8f0;\
            --shadow-sm: 0 1px 2px 0 rgba(0, 0, 0, 0.05);\
            --shadow-md: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);\
            --radius-md: 12px;\
            --radius-lg: 16px;\
        }\
        * { box-sizing: border-box; margin: 0; padding: 0; }\
        body {\
            font-family: 'Plus Jakarta Sans', 'Outfit', system-ui, sans-serif;\
            background-color: var(--bg-body);\
            color: var(--text-main);\
            display: flex;\
            min-height: 100vh;\
            overflow-x: hidden;\
        }\
        .sidebar {\
            width: 280px;\
            background: linear-gradient(180deg, var(--bg-sidebar) 0%, var(--bg-dark) 100%);\
            color: #fff;\
            padding: 32px 24px;\
            display: flex;\
            flex-direction: column;\
            flex-shrink: 0;\
            border-right: 1px solid rgba(255, 255, 255, 0.05);\
        }\
        .sidebar-brand {\
            display: flex;\
            align-items: center;\
            gap: 12px;\
            font-size: 1.45rem;\
            font-weight: 700;\
            letter-spacing: -0.025em;\
            margin-bottom: 36px;\
            color: #fff;\
            text-transform: uppercase;\
        }\
        .sidebar-brand span {\
            color: var(--primary);\
        }\
        .nav-menu {\
            list-style: none;\
            display: flex;\
            flex-direction: column;\
            gap: 8px;\
        }\
        .nav-item {\
            display: flex;\
            align-items: center;\
            gap: 12px;\
            padding: 12px 16px;\
            border-radius: var(--radius-md);\
            color: #94a3b8;\
            font-weight: 500;\
            text-decoration: none;\
            cursor: pointer;\
            transition: all 0.2s ease;\
        }\
        .nav-item:hover, .nav-item.active {\
            background-color: rgba(255, 255, 255, 0.06);\
            color: #fff;\
        }\
        .nav-item.active {\
            border-left: 4px solid var(--primary);\
            background-color: rgba(37, 99, 235, 0.15);\
        }\
        .main-content {\
            flex: 1;\
            padding: 40px;\
            overflow-y: auto;\
            display: flex;\
            flex-direction: column;\
            gap: 32px;\
        }\
        .header {\
            display: flex;\
            justify-content: space-between;\
            align-items: center;\
            border-bottom: 1px solid var(--border-color);\
            padding-bottom: 20px;\
        }\
        .header h1 {\
            font-size: 1.75rem;\
            font-weight: 700;\
            letter-spacing: -0.02em;\
            color: var(--bg-dark);\
        }\
        .header p {\
            color: var(--text-muted);\
            font-size: 0.9rem;\
            margin-top: 4px;\
        }\
        .sample-selector {\
            padding: 10px 16px;\
            font-family: inherit;\
            font-size: 0.95rem;\
            font-weight: 600;\
            color: var(--text-main);\
            background-color: #fff;\
            border: 1px solid var(--border-color);\
            border-radius: var(--radius-md);\
            outline: none;\
            cursor: pointer;\
            box-shadow: var(--shadow-sm);\
            transition: all 0.2s;\
        }\
        .sample-selector:focus {\
            border-color: var(--primary);\
        }\
        .tab-pane {\
            display: none;\
            flex-direction: column;\
            gap: 32px;\
        }\
        .tab-pane.active {\
            display: flex;\
        }\
        .kpi-grid {\
            display: grid;\
            grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));\
            gap: 20px;\
        }\
        .kpi-card {\
            background-color: var(--bg-card);\
            border: 1px solid var(--border-color);\
            border-radius: var(--radius-lg);\
            padding: 24px;\
            box-shadow: var(--shadow-sm);\
            display: flex;\
            flex-direction: column;\
            gap: 12px;\
            transition: transform 0.2s, box-shadow 0.2s;\
        }\
        .kpi-card:hover {\
            transform: translateY(-2px);\
            box-shadow: var(--shadow-md);\
        }\
        .kpi-title {\
            font-size: 0.85rem;\
            font-weight: 600;\
            text-transform: uppercase;\
            letter-spacing: 0.05em;\
            color: var(--text-muted);\
        }\
        .kpi-value {\
            font-size: 1.85rem;\
            font-weight: 700;\
            color: var(--bg-dark);\
        }\
        .kpi-subtitle {\
            font-size: 0.8rem;\
            color: var(--text-muted);\
        }\
        .card {\
            background-color: var(--bg-card);\
            border: 1px solid var(--border-color);\
            border-radius: var(--radius-lg);\
            padding: 28px;\
            box-shadow: var(--shadow-sm);\
            display: flex;\
            flex-direction: column;\
            gap: 20px;\
        }\
        .card-title {\
            font-size: 1.15rem;\
            font-weight: 600;\
            color: var(--bg-dark);\
            display: flex;\
            align-items: center;\
            justify-content: space-between;\
        }\
        .chart-container {\
            position: relative;\
            min-height: 320px;\
            width: 100%;\
        }\
        .chart-container canvas {\
            width: 100%;\
            height: 100%;\
            display: block;\
        }\
        .grid-2col {\
            display: grid;\
            grid-template-columns: repeat(auto-fit, minmax(480px, 1fr));\
            gap: 24px;\
        }\
        .table-container {\
            width: 100%;\
            overflow-x: auto;\
            border: 1px solid var(--border-color);\
            border-radius: var(--radius-md);\
        }\
        table {\
            width: 100%;\
            border-collapse: collapse;\
            font-size: 0.92rem;\
            text-align: left;\
        }\
        th {\
            background-color: #f8fafc;\
            padding: 14px 18px;\
            font-weight: 600;\
            color: var(--text-main);\
            border-bottom: 2px solid var(--border-color);\
        }\
        td {\
            padding: 14px 18px;\
            border-bottom: 1px solid var(--border-color);\
            color: var(--text-main);\
        }\
        tr:last-child td {\
            border-bottom: none;\
        }\
        .badge {\
            display: inline-flex;\
            align-items: center;\
            gap: 6px;\
            border-radius: 9999px;\
            background-color: var(--primary-light);\
            color: var(--primary);\
            padding: 4px 12px;\
            font-size: 0.78rem;\
            font-weight: 700;\
            text-transform: uppercase;\
            letter-spacing: 0.05em;\
        }\
        .search-bar {\
            padding: 12px 18px;\
            font-family: inherit;\
            font-size: 0.95rem;\
            border: 1px solid var(--border-color);\
            border-radius: var(--radius-md);\
            outline: none;\
            width: 100%;\
            max-width: 400px;\
            transition: all 0.2s;\
        }\
        .search-bar:focus {\
            border-color: var(--primary);\
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.15);\
        }\
        .links-flex {\
            display: flex;\
            gap: 12px;\
            flex-wrap: wrap;\
        }\
        .links-flex a {\
            display: inline-flex;\
            align-items: center;\
            gap: 8px;\
            text-decoration: none;\
            color: var(--primary);\
            background-color: var(--primary-light);\
            border: 1px solid var(--primary-border);\
            border-radius: var(--radius-md);\
            padding: 8px 16px;\
            font-size: 0.88rem;\
            font-weight: 600;\
            transition: all 0.2s;\
        }\
        .links-flex a:hover {\
            background-color: var(--primary);\
            color: #fff;\
        }\
        .snapshots-grid {\
            display: grid;\
            grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));\
            gap: 20px;\
        }\
        .snapshot-card {\
            border: 1px solid var(--border-color);\
            border-radius: var(--radius-md);\
            padding: 16px;\
            background-color: #f8fbff;\
            display: flex;\
            flex-direction: column;\
            gap: 10px;\
        }\
        .snapshot-card a {\
            color: var(--primary);\
            text-decoration: none;\
            font-weight: 700;\
            font-size: 0.95rem;\
        }\
        .snapshot-card .meta {\
            font-size: 0.82rem;\
            color: var(--text-muted);\
        }\
        pre {\
            background-color: var(--bg-dark);\
            color: #e2e8f0;\
            padding: 20px;\
            border-radius: var(--radius-md);\
            overflow-x: auto;\
            font-family: monospace;\
            font-size: 0.88rem;\
        }\
        </style>",
    );

    out.push_str("</head><body>");

    // Sidebar
    out.push_str("<div class=\"sidebar\">");
    out.push_str("<div class=\"sidebar-brand\">rs-qc<span>.engine</span></div>");
    out.push_str("<ul class=\"nav-menu\">");
    out.push_str("<li class=\"nav-item active\" data-tab=\"overview\" onclick=\"showTab('overview')\">Dashboard Overview</li>");
    out.push_str("<li class=\"nav-item\" id=\"nav-fastq\" data-tab=\"fastq\" onclick=\"showTab('fastq')\" style=\"display:none;\">FASTQ Metrics</li>");
    out.push_str("<li class=\"nav-item\" id=\"nav-align\" data-tab=\"align\" onclick=\"showTab('align')\" style=\"display:none;\">Alignment Metrics</li>");
    out.push_str("<li class=\"nav-item\" id=\"nav-dna\" data-tab=\"dna\" onclick=\"showTab('dna')\" style=\"display:none;\">DNA Coverage</li>");
    out.push_str("<li class=\"nav-item\" id=\"nav-rna\" data-tab=\"rna\" onclick=\"showTab('rna')\" style=\"display:none;\">RNA Read Dist</li>");
    out.push_str("<li class=\"nav-item\" data-tab=\"raw-data\" onclick=\"showTab('raw-data')\">Search & Flat Table</li>");
    out.push_str("</ul>");
    out.push_str("</div>");

    // Main Pane
    out.push_str("<div class=\"main-content\">");
    
    // Header
    out.push_str("<div class=\"header\">");
    out.push_str("<div>");
    out.push_str(&format!("<h1>{}</h1>", escape_html(&document.title)));
    out.push_str(&format!("<p>Report generated by rs-qc v{}</p>", escape_html(&document.generated_by)));
    out.push_str("</div>");
    out.push_str("<select class=\"sample-selector\" id=\"sampleSelect\" onchange=\"onSampleChange()\"></select>");
    out.push_str("</div>");

    // Tabs Container
    // Tab: Overview
    out.push_str("<div class=\"tab-pane active\" id=\"tab-overview\">");
    out.push_str("<div class=\"kpi-grid\" id=\"overviewKpis\"></div>");
    out.push_str("<div class=\"grid-2col\">");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Active QC Modules</div><div class=\"table-container\"><table>");
    out.push_str("<thead><tr><th>Module</th><th>Sample</th><th>Source JSON File</th></tr></thead>");
    out.push_str("<tbody id=\"overviewModules\"></tbody></table></div></div>");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Outputs & Export Files</div><div class=\"links-flex\" id=\"overviewLinks\"></div></div>");
    out.push_str("</div></div>");

    // Tab: FASTQ
    out.push_str("<div class=\"tab-pane\" id=\"tab-fastq\">");
    out.push_str("<div class=\"grid-2col\">");
    out.push_str("<div class=\"card\"><div class=\"card-title\">GC Content Distribution</div><div class=\"chart-container\"><canvas id=\"fastqGcChart\"></canvas></div></div>");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Mean Quality Distribution</div><div class=\"chart-container\"><canvas id=\"fastqQualChart\"></canvas></div></div>");
    out.push_str("</div></div>");

    // Tab: Alignment
    out.push_str("<div class=\"tab-pane\" id=\"tab-align\">");
    out.push_str("<div class=\"grid-2col\">");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Insert Size Distribution</div><div class=\"chart-container\"><canvas id=\"alignInsertChart\"></canvas></div></div>");
    out.push_str("<div class=\"card\"><div class=\"card-title\">MAPQ Score Distribution</div><div class=\"chart-container\"><canvas id=\"alignMapqChart\"></canvas></div></div>");
    out.push_str("</div></div>");

    // Tab: DNA
    out.push_str("<div class=\"tab-pane\" id=\"tab-dna\">");
    out.push_str("<div class=\"grid-2col\">");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Depth Histogram</div><div class=\"chart-container\"><canvas id=\"dnaDepthChart\"></canvas></div></div>");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Breadth of Coverage</div><div class=\"table-container\"><table><thead><tr><th>Breadth threshold</th><th>Coverage Percentage</th></tr></thead><tbody id=\"dnaBreadthTable\"></tbody></table></div></div>");
    out.push_str("</div></div>");

    // Tab: RNA
    out.push_str("<div class=\"tab-pane\" id=\"tab-rna\">");
    out.push_str("<div class=\"grid-2col\">");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Read Distribution Profile</div><div class=\"chart-container\"><canvas id=\"rnaDistChart\"></canvas></div></div>");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Junction Annotation Summary</div><div class=\"table-container\"><table><thead><tr><th>Junction Class</th><th>Unique Count</th><th>Total Reads</th></tr></thead><tbody id=\"rnaJunctionTable\"></tbody></table></div></div>");
    out.push_str("</div>");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Clinical Genomic Snapshots</div><div class=\"snapshots-grid\" id=\"rnaSnapshotsGrid\"></div></div>");
    out.push_str("</div>");

    // Tab: Raw Data
    out.push_str("<div class=\"tab-pane\" id=\"tab-raw-data\">");
    out.push_str("<div class=\"card\">");
    out.push_str("<div class=\"card-title\">Instant QC Metric Search & Filter <input class=\"search-bar\" id=\"searchBar\" placeholder=\"Search metric key or value...\" oninput=\"onSearch()\"></input></div>");
    out.push_str("<div class=\"table-container\">");
    out.push_str("<table><thead><tr><th>Metric Key</th><th>Metric Value</th></tr></thead><tbody id=\"flatMetricsTable\"></tbody></table>");
    out.push_str("</div></div>");
    out.push_str("<div class=\"card\"><div class=\"card-title\">Raw Summary Metrics JSON Document</div><pre id=\"rawJsonDisplay\"></pre></div>");
    out.push_str("</div>");

    out.push_str("</div>"); // Close main-content

    // Inject raw JSON data
    out.push_str(&format!(
        "<script>const reportData = {};</script>",
        safe_json_data
    ));

    // Dashboard controller logic
    out.push_str(
        "<script>\
        let charts = {};\
        \
        function initDashboard() {\
            const select = document.getElementById('sampleSelect');\
            const samples = [...new Set(reportData.sections.map(s => s.sample))];\
            samples.forEach(s => {\
                const opt = document.createElement('option');\
                opt.value = s;\
                opt.textContent = s;\
                select.appendChild(opt);\
            });\
            \
            // Check which modules are present overall\
            const modules = new Set(reportData.sections.map(s => s.module));\
            if (modules.has('fastq')) document.getElementById('nav-fastq').style.display = 'flex';\
            if (modules.has('align')) document.getElementById('nav-align').style.display = 'flex';\
            if (modules.has('dna')) document.getElementById('nav-dna').style.display = 'flex';\
            if (modules.has('rna')) document.getElementById('nav-rna').style.display = 'flex';\
            \
            onSampleChange();\
        }\
        \
        function showTab(tabId) {\
            document.querySelectorAll('.nav-item').forEach(item => item.classList.remove('active'));\
            document.querySelectorAll('.tab-pane').forEach(pane => pane.classList.remove('active'));\
            \
            const activeNav = document.querySelector(`.nav-item[data-tab=\"${tabId}\"]`);\
            if (activeNav) activeNav.classList.add('active');\
            \
            const activePane = document.getElementById('tab-' + tabId);\
            if (activePane) activePane.classList.add('active');\
            \
            // Re-render active tab charts if needed\
            renderActiveTabCharts(tabId);\
        }\
        \
        function onSampleChange() {\
            const sample = document.getElementById('sampleSelect').value;\
            const sections = reportData.sections.filter(s => s.sample === sample);\
            \
            // 1. Populate Overview Modules Table\
            const modTbody = document.getElementById('overviewModules');\
            modTbody.innerHTML = '';\
            sections.forEach(s => {\
                const tr = document.createElement('tr');\
                tr.innerHTML = `<td><span class=\"badge\">${s.module}</span></td><td>${s.sample}</td><td>${s.source}</td>`;\
                modTbody.appendChild(tr);\
            });\
            \
            // 2. Populate Overview Links (dummy files based on prefixes)\
            const linksDiv = document.getElementById('overviewLinks');\
            linksDiv.innerHTML = '';\
            sections.forEach(s => {\
                const baseName = s.source.substring(s.source.lastIndexOf('/') + 1);\
                const prefix = baseName.split('.')[0];\
                const files = [\
                    { label: 'Summary JSON', ext: '.summary.json' },\
                    { label: 'Summary TXT/TSV', ext: s.module === 'dna' ? '.dna.summary.tsv' : s.module === 'fastq' ? '.fastq.summary.txt' : s.module === 'align' ? '.align.summary.txt' : '.rna_qc.txt' }\
                ];\
                files.forEach(f => {\
                    const a = document.createElement('a');\
                    a.href = `./${prefix}${f.ext}`;\
                    a.target = '_blank';\
                    a.textContent = `${s.module.toUpperCase()} - ${f.label}`;\
                    linksDiv.appendChild(a);\
                });\
            });\
            \
            // 3. Build KPI Cards\
            buildKpis(sections);\
            \
            // 4. Reset Charts\
            Object.values(charts).forEach(c => c.destroy());\
            charts = {};\
            \
            // 5. Populate Detailed Raw Table & JSON\
            populateFlatMetrics(sections);\
            document.getElementById('rawJsonDisplay').textContent = JSON.stringify(sections.map(s => s.metrics), null, 2);\
            \
            // 6. Draw current active tab charts\
            const activeTab = document.querySelector('.nav-item.active').getAttribute('data-tab');\
            renderActiveTabCharts(activeTab);\
        }\
        \
        function buildKpis(sections) {\
            const grid = document.getElementById('overviewKpis');\
            grid.innerHTML = '';\
            \
            sections.forEach(s => {\
                const m = s.metrics;\
                if (s.module === 'fastq') {\
                    addKpi(grid, 'Total Reads', formatNum(m.total_reads), 'FASTQ QC');\
                    addKpi(grid, 'Total Bases', formatNum(m.total_bases), 'FASTQ QC');\
                    if (m.duplication_estimate) {\
                        addKpi(grid, 'Est. Duplication Rate', (m.duplication_estimate * 100).toFixed(1) + '%', 'FASTQ QC');\
                    }\
                } else if (s.module === 'align') {\
                    addKpi(grid, 'Mapped Records', formatNum(m.mapped_records) + ` (${(m.mapped_records/m.total_records*100).toFixed(1)}%)`, 'Alignment QC');\
                    if (m.properly_paired_records) {\
                        addKpi(grid, 'Properly Paired', formatNum(m.properly_paired_records), 'Alignment QC');\
                    }\
                    if (m.accuracy && m.accuracy.accuracy_sum) {\
                        const acc = (m.accuracy.accuracy_sum / m.accuracy.records_with_de * 100).toFixed(2);\
                        addKpi(grid, 'Alignment Accuracy', acc + '%', 'Alignment QC');\
                    }\
                } else if (s.module === 'dna') {\
                    addKpi(grid, 'Mean Target Depth', m.mean_depth.toFixed(1) + 'x', 'DNA Coverage');\
                    addKpi(grid, 'Breadth at 10x', (m.breadth_10x * 100).toFixed(1) + '%', 'DNA Coverage');\
                    addKpi(grid, 'Breadth at 30x', (m.breadth_30x * 100).toFixed(1) + '%', 'DNA Coverage');\
                } else if (s.module === 'rna') {\
                    addKpi(grid, 'Aligned QC Reads', formatNum(m.aligned_qc_reads), 'RNA QC');\
                    addKpi(grid, 'Strandness', m.inferred_strandness, 'RNA QC');\
                    addKpi(grid, 'Active Genes', formatNum(m.active_genes), 'RNA QC');\
                }\
            });\
        }\
        \
        function addKpi(container, title, value, module) {\
            const card = document.createElement('div');\
            card.className = 'kpi-card';\
            card.innerHTML = `<span class=\"kpi-title\">${title}</span><span class=\"kpi-value\">${value}</span><span class=\"kpi-subtitle\">${module}</span>`;\
            container.appendChild(card);\
        }\
        \
        function formatNum(n) {\
            if (n === undefined || n === null) return 'NA';\
            return Number(n).toLocaleString();\
        }\
        \
        function populateFlatMetrics(sections) {\
            const tbody = document.getElementById('flatMetricsTable');\
            tbody.innerHTML = '';\
            sections.forEach(s => {\
                const flat = flattenObj(s.metrics, s.module);\
                Object.entries(flat).forEach(([k, v]) => {\
                    const tr = document.createElement('tr');\
                    tr.innerHTML = `<td class=\"key\" style=\"font-weight:600;\">${k}</td><td>${v}</td>`;\
                    tbody.appendChild(tr);\
                });\
            });\
        }\
        \
        function flattenObj(val, prefix) {\
            let res = {};\
            if (typeof val === 'object' && val !== null) {\
                if (Array.isArray(val)) {\
                    res[prefix] = `${val.length} items`;\
                } else {\
                    Object.entries(val).forEach(([k, child]) => {\
                        if (k === 'length_hist' || k === 'gc_hist' || k === 'mean_quality_hist' || k === 'mapq_hist' || k === 'read_length_hist' || k === 'insert_size_hist' || k === 'depth_hist' || k === 'accuracy_hist') {\
                            // Skip huge histograms in detailed flat table view\
                            return;\
                        }\
                        Object.assign(res, flattenObj(child, `${prefix}.${k}`));\
                    });\
                }\
            } else {\
                res[prefix] = val === null ? 'NA' : val;\
            }\
            return res;\
        }\
        \
        function onSearch() {\
            const query = document.getElementById('searchBar').value.toLowerCase();\
            document.querySelectorAll('#flatMetricsTable tr').forEach(tr => {\
                const text = tr.textContent.toLowerCase();\
                tr.style.display = text.includes(query) ? '' : 'none';\
            });\
        }\
        \
        function renderActiveTabCharts(tabId) {
            if (typeof Chart === 'undefined') {
                console.warn('Chart.js is not loaded (running offline). Charts will not be drawn.');
                return;
            }
            const sample = document.getElementById('sampleSelect').value;
            const sections = reportData.sections.filter(s => s.sample === sample);\
            \
            if (tabId === 'fastq') {\
                const sec = sections.find(s => s.module === 'fastq');\
                if (sec) drawFastqCharts(sec.metrics);\
            } else if (tabId === 'align') {\
                const sec = sections.find(s => s.module === 'align');\
                if (sec) drawAlignCharts(sec.metrics);\
            } else if (tabId === 'dna') {\
                const sec = sections.find(s => s.module === 'dna');\
                if (sec) drawDnaCharts(sec.metrics);\
            } else if (tabId === 'rna') {\
                const sec = sections.find(s => s.module === 'rna');\
                if (sec) drawRnaCharts(sec.metrics);\
            }\
        }\
        \
        function drawFastqCharts(m) {\
            if (m.gc_hist && !charts['fastqGc']) {\
                const ctx = document.getElementById('fastqGcChart').getContext('2d');\
                const labels = Object.keys(m.gc_hist);\
                const values = Object.values(m.gc_hist);\
                charts['fastqGc'] = new Chart(ctx, {\
                    type: 'line',\
                    data: {\
                        labels,\
                        datasets: [{\
                            label: 'GC Percentage Count',\
                            data: values,\
                            borderColor: '#2563eb',\
                            backgroundColor: 'rgba(37,99,235,0.05)',\
                            borderWidth: 2,\
                            fill: true,\
                            tension: 0.3\
                        }]\
                    },\
                    options: { responsive: true, maintainAspectRatio: false }\
                });\
            }\
            if (m.mean_quality_hist && !charts['fastqQual']) {\
                const ctx = document.getElementById('fastqQualChart').getContext('2d');\
                const labels = Object.keys(m.mean_quality_hist);\
                const values = Object.values(m.mean_quality_hist);\
                charts['fastqQual'] = new Chart(ctx, {\
                    type: 'bar',\
                    data: {\
                        labels,\
                        datasets: [{\
                            label: 'Mean Quality Score Count',\
                            data: values,\
                            backgroundColor: '#10b981'\
                        }]\
                    },\
                    options: { responsive: true, maintainAspectRatio: false }\
                });\
            }\
        }\
        \
        function drawAlignCharts(m) {\
            if (m.insert_size_hist && !charts['alignInsert']) {\
                const ctx = document.getElementById('alignInsertChart').getContext('2d');\
                const labels = Object.keys(m.insert_size_hist);\
                const values = Object.values(m.insert_size_hist);\
                charts['alignInsert'] = new Chart(ctx, {\
                    type: 'line',\
                    data: {\
                        labels,\
                        datasets: [{\
                            label: 'Insert Size Count',\
                            data: values,\
                            borderColor: '#2563eb',\
                            backgroundColor: 'rgba(37,99,235,0.05)',\
                            borderWidth: 2,\
                            fill: true,\
                            tension: 0.3\
                        }]\
                    },\
                    options: { responsive: true, maintainAspectRatio: false }\
                });\
            }\
            if (m.mapq_hist && !charts['alignMapq']) {\
                const ctx = document.getElementById('alignMapqChart').getContext('2d');\
                const labels = Object.keys(m.mapq_hist);\
                const values = Object.values(m.mapq_hist);\
                charts['alignMapq'] = new Chart(ctx, {\
                    type: 'bar',\
                    data: {\
                        labels,\
                        datasets: [{\
                            label: 'MAPQ Count',\
                            data: values,\
                            backgroundColor: '#1e293b'\
                        }]\
                    },\
                    options: { responsive: true, maintainAspectRatio: false }\
                });\
            }\
        }\
        \
        function drawDnaCharts(m) {\
            if (m.depth_hist && !charts['dnaDepth']) {\
                const ctx = document.getElementById('dnaDepthChart').getContext('2d');\
                const labels = Object.keys(m.depth_hist);\
                const values = Object.values(m.depth_hist);\
                charts['dnaDepth'] = new Chart(ctx, {\
                    type: 'bar',\
                    data: {\
                        labels,\
                        datasets: [{\
                            label: 'Depth Count',\
                            data: values,\
                            backgroundColor: '#2563eb'\
                        }]\
                    },\
                    options: { responsive: true, maintainAspectRatio: false }\
                });\
            }\
            \
            // Populate DNA Breadth table\
            const tbody = document.getElementById('dnaBreadthTable');\
            tbody.innerHTML = `\
                <tr><td>1x</td><td>${(m.breadth_1x*100).toFixed(2)}%</td></tr>\
                <tr><td>5x</td><td>${(m.breadth_5x*100).toFixed(2)}%</td></tr>\
                <tr><td>10x</td><td>${(m.breadth_10x*100).toFixed(2)}%</td></tr>\
                <tr><td>20x</td><td>${(m.breadth_20x*100).toFixed(2)}%</td></tr>\
                <tr><td>30x</td><td>${(m.breadth_30x*100).toFixed(2)}%</td></tr>\
            `;\
        }\
        \
        function drawRnaCharts(m) {\
            if (!charts['rnaDist']) {\
                const ctx = document.getElementById('rnaDistChart').getContext('2d');\
                charts['rnaDist'] = new Chart(ctx, {\
                    type: 'doughnut',\
                    data: {\
                        labels: ['Exonic', 'Intronic', 'Flanking', 'Intergenic'],\
                        datasets: [{\
                            data: [m.exonic_reads, m.intronic_reads, m.flanking_reads, m.intergenic_reads],\
                            backgroundColor: ['#10b981', '#f59e0b', '#3b82f6', '#94a3b8']\
                        }]\
                    },\
                    options: { responsive: true, maintainAspectRatio: false }\
                });\
            }\
            \
            // Populate Junction Table\
            const tbody = document.getElementById('rnaJunctionTable');\
            tbody.innerHTML = `\
                <tr><td>Known Junctions</td><td>${m.unique_known_junctions || 0}</td><td>${m.total_known_junction_reads || 0}</td></tr>\
                <tr><td>Partial Novel Junctions</td><td>${m.unique_partial_novel_junctions || 0}</td><td>${m.total_partial_novel_junction_reads || 0}</td></tr>\
                <tr><td>Fully Novel Junctions</td><td>${m.unique_novel_junctions || 0}</td><td>${m.total_novel_junction_reads || 0}</td></tr>\
            `;\
            \
            // Render Clinical Genomic Snapshots gallery if available\
            const snapGrid = document.getElementById('rnaSnapshotsGrid');\
            snapGrid.innerHTML = '';\
            const baseName = reportData.sections.find(s => s.module === 'rna').source.substring(reportData.sections.find(s => s.module === 'rna').source.lastIndexOf('/') + 1);\
            const prefix = baseName.split('.')[0];\
            const snaps = [\
                { gene: 'GAPDH', file: `${prefix}.GAPDH.png`, reason: 'Housekeeping Gene Body' }\
            ];\
            snaps.forEach(s => {\
                const card = document.createElement('div');\
                card.className = 'snapshot-card';\
                card.innerHTML = `<a href=\"./snapshots/${s.file}\" target=\"_blank\">${s.gene} Snapshot</a><span class=\"meta\">${s.reason}</span>`;\
                snapGrid.appendChild(card);\
            });\
        }\
        \
        window.onload = initDashboard;\
        </script>",
    );

    out.push_str(
        r##"<script>
        function fitCanvas(canvas) {
            const rect = canvas.getBoundingClientRect();
            const width = Math.max(1, Math.round(rect.width || canvas.clientWidth || 300));
            const height = Math.max(1, Math.round(rect.height || canvas.clientHeight || 240));
            const dpr = window.devicePixelRatio || 1;
            canvas.width = Math.round(width * dpr);
            canvas.height = Math.round(height * dpr);
            const ctx = canvas.getContext('2d');
            ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            return { ctx, width, height };
        }

        function clearCanvas(canvas) {
            const { ctx, width, height } = fitCanvas(canvas);
            ctx.clearRect(0, 0, width, height);
            return { ctx, width, height };
        }

        function drawGridAxes(ctx, width, height, margin, maxValue, tickCount) {
            ctx.save();
            ctx.strokeStyle = '#e2e8f0';
            ctx.fillStyle = '#64748b';
            ctx.lineWidth = 1;
            ctx.font = '12px sans-serif';
            ctx.textBaseline = 'middle';

            const plotWidth = width - margin.left - margin.right;
            const plotHeight = height - margin.top - margin.bottom;
            const safeMax = Math.max(maxValue, 1);

            for (let i = 0; i <= tickCount; i++) {
                const ratio = i / tickCount;
                const y = margin.top + plotHeight - ratio * plotHeight;
                ctx.beginPath();
                ctx.moveTo(margin.left, y);
                ctx.lineTo(width - margin.right, y);
                ctx.stroke();
                ctx.fillText(Math.round(safeMax * ratio).toLocaleString(), 8, y);
            }

            ctx.strokeStyle = '#0f172a';
            ctx.beginPath();
            ctx.moveTo(margin.left, margin.top);
            ctx.lineTo(margin.left, height - margin.bottom);
            ctx.lineTo(width - margin.right, height - margin.bottom);
            ctx.stroke();

            ctx.restore();
        }

        function drawLineFallback(canvasId, labels, values, color, fill) {
            const canvas = document.getElementById(canvasId);
            if (!canvas) return;
            const { ctx, width, height } = clearCanvas(canvas);
            const margin = { top: 18, right: 18, bottom: 36, left: 48 };
            const plotWidth = width - margin.left - margin.right;
            const plotHeight = height - margin.top - margin.bottom;
            const nums = values.map(v => Number(v) || 0);
            const maxValue = Math.max(1, ...nums);
            const step = nums.length > 12 ? Math.ceil(nums.length / 10) : 1;

            drawGridAxes(ctx, width, height, margin, maxValue, 5);

            if (nums.length === 0) return;

            const points = nums.map((value, index) => {
                const x = margin.left + (nums.length === 1 ? plotWidth / 2 : (index / (nums.length - 1)) * plotWidth);
                const y = margin.top + plotHeight - (value / maxValue) * plotHeight;
                return { x, y, value };
            });

            if (fill) {
                ctx.beginPath();
                ctx.moveTo(points[0].x, margin.top + plotHeight);
                points.forEach((point, index) => {
                    if (index === 0) {
                        ctx.lineTo(point.x, point.y);
                    } else {
                        ctx.lineTo(point.x, point.y);
                    }
                });
                ctx.lineTo(points[points.length - 1].x, margin.top + plotHeight);
                ctx.closePath();
                ctx.fillStyle = 'rgba(37, 99, 235, 0.08)';
                ctx.fill();
            }

            ctx.beginPath();
            ctx.strokeStyle = color;
            ctx.lineWidth = 2;
            points.forEach((point, index) => {
                if (index === 0) {
                    ctx.moveTo(point.x, point.y);
                } else {
                    ctx.lineTo(point.x, point.y);
                }
            });
            ctx.stroke();

            ctx.fillStyle = color;
            points.forEach(point => {
                ctx.beginPath();
                ctx.arc(point.x, point.y, 1.8, 0, Math.PI * 2);
                ctx.fill();
            });

            ctx.save();
            ctx.fillStyle = '#64748b';
            ctx.font = '11px sans-serif';
            ctx.textBaseline = 'top';
            ctx.textAlign = 'center';
            labels.forEach((label, index) => {
                if (index % step !== 0 && index !== labels.length - 1) return;
                const x = points[index] ? points[index].x : margin.left;
                ctx.fillText(String(label), x, height - margin.bottom + 8);
            });
            ctx.restore();
        }

        function drawBarFallback(canvasId, labels, values, color) {
            const canvas = document.getElementById(canvasId);
            if (!canvas) return;
            const { ctx, width, height } = clearCanvas(canvas);
            const margin = { top: 18, right: 18, bottom: 40, left: 48 };
            const plotWidth = width - margin.left - margin.right;
            const plotHeight = height - margin.top - margin.bottom;
            const nums = values.map(v => Number(v) || 0);
            const maxValue = Math.max(1, ...nums);
            const step = nums.length > 12 ? Math.ceil(nums.length / 10) : 1;
            const barWidth = nums.length > 0 ? plotWidth / nums.length : plotWidth;

            drawGridAxes(ctx, width, height, margin, maxValue, 5);

            ctx.fillStyle = color;
            nums.forEach((value, index) => {
                const scaledHeight = (value / maxValue) * plotHeight;
                const x = margin.left + index * barWidth + Math.max(0, barWidth * 0.12);
                const y = margin.top + plotHeight - scaledHeight;
                ctx.fillRect(x, y, Math.max(1, barWidth * 0.76), scaledHeight);
            });

            ctx.save();
            ctx.fillStyle = '#64748b';
            ctx.font = '11px sans-serif';
            ctx.textBaseline = 'top';
            ctx.textAlign = 'center';
            labels.forEach((label, index) => {
                if (index % step !== 0 && index !== labels.length - 1) return;
                const x = margin.left + index * barWidth + barWidth / 2;
                ctx.fillText(String(label), x, height - margin.bottom + 8);
            });
            ctx.restore();
        }

        function drawDoughnutFallback(canvasId, labels, values, colors) {
            const canvas = document.getElementById(canvasId);
            if (!canvas) return;
            const { ctx, width, height } = clearCanvas(canvas);
            const nums = values.map(v => Math.max(0, Number(v) || 0));
            const total = nums.reduce((sum, value) => sum + value, 0) || 1;
            const radius = Math.max(36, Math.min(width, height) * 0.22);
            const innerRadius = radius * 0.58;
            const centerX = Math.min(width * 0.38, width / 2 - 24);
            const centerY = height / 2;
            let start = -Math.PI / 2;

            nums.forEach((value, index) => {
                const slice = (value / total) * Math.PI * 2;
                ctx.beginPath();
                ctx.moveTo(centerX, centerY);
                ctx.arc(centerX, centerY, radius, start, start + slice);
                ctx.closePath();
                ctx.fillStyle = colors[index % colors.length];
                ctx.fill();
                start += slice;
            });

            ctx.globalCompositeOperation = 'destination-out';
            ctx.beginPath();
            ctx.arc(centerX, centerY, innerRadius, 0, Math.PI * 2);
            ctx.fill();
            ctx.globalCompositeOperation = 'source-over';

            const legendX = Math.max(centerX + radius + 18, width * 0.58);
            const legendY = Math.max(18, height / 2 - (labels.length * 20) / 2);
            ctx.font = '12px sans-serif';
            ctx.textBaseline = 'middle';
            labels.forEach((label, index) => {
                const y = legendY + index * 22;
                ctx.fillStyle = colors[index % colors.length];
                ctx.fillRect(legendX, y - 6, 12, 12);
                ctx.fillStyle = '#1e293b';
                ctx.fillText(String(label), legendX + 18, y);
            });
        }

        function drawFastqChartsFallback(m) {
            if (m.gc_hist) {
                drawLineFallback('fastqGcChart', Object.keys(m.gc_hist), Object.values(m.gc_hist), '#2563eb', true);
            }
            if (m.mean_quality_hist) {
                drawBarFallback('fastqQualChart', Object.keys(m.mean_quality_hist), Object.values(m.mean_quality_hist), '#10b981');
            }
        }

        function drawAlignChartsFallback(m) {
            if (m.insert_size_hist) {
                drawLineFallback('alignInsertChart', Object.keys(m.insert_size_hist), Object.values(m.insert_size_hist), '#2563eb', true);
            }
            if (m.mapq_hist) {
                drawBarFallback('alignMapqChart', Object.keys(m.mapq_hist), Object.values(m.mapq_hist), '#1e293b');
            }
        }

        function drawDnaChartsFallback(m) {
            if (m.depth_hist) {
                drawBarFallback('dnaDepthChart', Object.keys(m.depth_hist), Object.values(m.depth_hist), '#2563eb');
            }
            const tbody = document.getElementById('dnaBreadthTable');
            tbody.innerHTML = `
                <tr><td>1x</td><td>${(m.breadth_1x * 100).toFixed(2)}%</td></tr>
                <tr><td>5x</td><td>${(m.breadth_5x * 100).toFixed(2)}%</td></tr>
                <tr><td>10x</td><td>${(m.breadth_10x * 100).toFixed(2)}%</td></tr>
                <tr><td>20x</td><td>${(m.breadth_20x * 100).toFixed(2)}%</td></tr>
                <tr><td>30x</td><td>${(m.breadth_30x * 100).toFixed(2)}%</td></tr>
            `;
        }

        function drawRnaChartsFallback(m) {
            drawDoughnutFallback(
                'rnaDistChart',
                ['Exonic', 'Intronic', 'Flanking', 'Intergenic'],
                [m.exonic_reads, m.intronic_reads, m.flanking_reads, m.intergenic_reads],
                ['#10b981', '#f59e0b', '#3b82f6', '#94a3b8']
            );

            const tbody = document.getElementById('rnaJunctionTable');
            tbody.innerHTML = `
                <tr><td>Known Junctions</td><td>${m.unique_known_junctions || 0}</td><td>${m.total_known_junction_reads || 0}</td></tr>
                <tr><td>Partial Novel Junctions</td><td>${m.unique_partial_novel_junctions || 0}</td><td>${m.total_partial_novel_junction_reads || 0}</td></tr>
                <tr><td>Fully Novel Junctions</td><td>${m.unique_novel_junctions || 0}</td><td>${m.total_novel_junction_reads || 0}</td></tr>
            `;

            const snapGrid = document.getElementById('rnaSnapshotsGrid');
            snapGrid.innerHTML = '';
            const rnaSection = reportData.sections.find(s => s.module === 'rna');
            if (!rnaSection) return;
            const baseName = rnaSection.source.substring(rnaSection.source.lastIndexOf('/') + 1);
            const prefix = baseName.split('.')[0];
            const snaps = [
                { gene: 'GAPDH', file: `${prefix}.GAPDH.png`, reason: 'Housekeeping Gene Body' }
            ];
            snaps.forEach(s => {
                const card = document.createElement('div');
                card.className = 'snapshot-card';
                card.innerHTML = `<a href="./snapshots/${s.file}" target="_blank">${s.gene} Snapshot</a><span class="meta">${s.reason}</span>`;
                snapGrid.appendChild(card);
            });
        }

        function renderActiveTabCharts(tabId) {
            const sample = document.getElementById('sampleSelect').value;
            const sections = reportData.sections.filter(s => s.sample === sample);

            if (typeof Chart === 'undefined') {
                if (tabId === 'fastq') {
                    const sec = sections.find(s => s.module === 'fastq');
                    if (sec) drawFastqChartsFallback(sec.metrics);
                } else if (tabId === 'align') {
                    const sec = sections.find(s => s.module === 'align');
                    if (sec) drawAlignChartsFallback(sec.metrics);
                } else if (tabId === 'dna') {
                    const sec = sections.find(s => s.module === 'dna');
                    if (sec) drawDnaChartsFallback(sec.metrics);
                } else if (tabId === 'rna') {
                    const sec = sections.find(s => s.module === 'rna');
                    if (sec) drawRnaChartsFallback(sec.metrics);
                }
                return;
            }

            if (tabId === 'fastq') {
                const sec = sections.find(s => s.module === 'fastq');
                if (sec) drawFastqCharts(sec.metrics);
            } else if (tabId === 'align') {
                const sec = sections.find(s => s.module === 'align');
                if (sec) drawAlignCharts(sec.metrics);
            } else if (tabId === 'dna') {
                const sec = sections.find(s => s.module === 'dna');
                if (sec) drawDnaCharts(sec.metrics);
            } else if (tabId === 'rna') {
                const sec = sections.find(s => s.module === 'rna');
                if (sec) drawRnaCharts(sec.metrics);
            }
        }
        </script>"##,
    );

    out.push_str("</body></html>");
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
    let Some(prefix) = summary_prefix_from_path(source) else {
        return Vec::new();
    };
    let dir = source.parent().unwrap_or_else(|| Path::new("."));
    let candidates = match section.module.as_str() {
        "fastq" => vec![
            format!("{prefix}.fastq.summary.json"),
            format!("{prefix}.fastq.summary.txt"),
            format!("{prefix}.fastq.per_base.tsv"),
            format!("{prefix}.fastq.length_distribution.tsv"),
            format!("{prefix}.fastq.gc_distribution.tsv"),
            format!("{prefix}.fastq.mean_quality_distribution.tsv"),
            format!("{prefix}.fastq.overrepresented.tsv"),
            format!("{prefix}.fastq.kmers.tsv"),
        ],
        "align" => vec![
            format!("{prefix}.align.summary.json"),
            format!("{prefix}.align.summary.txt"),
            format!("{prefix}.align.mapq.tsv"),
            format!("{prefix}.align.read_length.tsv"),
            format!("{prefix}.align.insert_size.tsv"),
            format!("{prefix}.align.cigar.tsv"),
            format!("{prefix}.align.contigs.tsv"),
            format!("{prefix}.align.de_accuracy.tsv"),
        ],
        "dna" => vec![
            format!("{prefix}.dna.summary.json"),
            format!("{prefix}.dna.summary.tsv"),
            format!("{prefix}.dna.depth_hist.tsv"),
            format!("{prefix}.dna.contigs.tsv"),
            format!("{prefix}.dna.windows.tsv"),
            format!("{prefix}.dna.targets.tsv"),
        ],
        "rna" => vec![
            format!("{prefix}.rna.summary.json"),
            format!("{prefix}.rna_qc.txt"),
            format!("{prefix}.inner_distance.tsv"),
            format!("{prefix}.geneBodyCoverage.txt"),
            format!("{prefix}.geneBodyCoverage.svg"),
            format!("{prefix}.rna.qc_summary.svg"),
            format!("{prefix}.rna_snapshots.tsv"),
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

fn rna_summary_figure_link(section: &ReportSection) -> Option<(String, String)> {
    let source = Path::new(&section.source);
    let prefix = summary_prefix_from_path(source)?;
    let dir = source.parent().unwrap_or_else(|| Path::new("."));
    let name = format!("{prefix}.{}.rna.qc_summary.svg", section.sample);
    let path = dir.join(&name);
    if path.exists() {
        Some((name, format!("RNA QC summary for {}", section.sample)))
    } else {
        None
    }
}

fn rna_snapshot_items(section: &ReportSection) -> Option<Vec<(String, String, String)>> {
    let source = Path::new(&section.source);
    let prefix = summary_prefix_from_path(source)?;
    let dir = source.parent().unwrap_or_else(|| Path::new("."));
    let manifest = dir.join(format!("{prefix}.{}.rna_snapshots.tsv", section.sample));
    let raw = fs::read_to_string(manifest).ok()?;
    let mut items = Vec::new();
    for line in raw.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 7 {
            continue;
        }
        let gene = cols[0].to_string();
        let snapshot_path = cols[5].to_string();
        let reason = cols[6].to_string();
        if snapshot_path.is_empty() {
            continue;
        }
        items.push((gene, snapshot_path, reason));
    }
    Some(items)
}

fn summary_prefix_from_path(path: &Path) -> Option<String> {
    let name = path.file_name()?.to_str()?;
    for suffix in [
        ".fastq.summary.json",
        ".align.summary.json",
        ".rna.summary.json",
        ".dna.summary.json",
        ".summary.json",
    ] {
        if let Some(prefix) = name.strip_suffix(suffix) {
            return Some(prefix.to_string());
        }
    }
    None
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
    } else if name.contains(".dna.summary.json") {
        "dna"
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
        let dna_path = tmp.path().join("sample.dna.summary.json");
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
            dna_path.to_str().unwrap(),
            "dna",
            "sample",
            &json!({"total_records": 40, "mean_depth": 30.5}),
        )
        .expect("write dna summary");
        write_summary_json(
            rna_path.to_str().unwrap(),
            "rna",
            "sample",
            &json!({"aligned_qc_reads": 5, "mtdna_reads": 1}),
        )
        .expect("write rna summary");

        let files = resolve_input_files(&["sample".to_string()]).expect("resolve prefix");
        assert_eq!(files.len(), 4);

        let document = build_document(&files).expect("build document");
        assert_eq!(document.sections.len(), 4);
        assert!(document.sections.iter().any(|s| s.module == "fastq"));
        assert!(document.sections.iter().any(|s| s.module == "align"));
        assert!(document.sections.iter().any(|s| s.module == "dna"));
        assert!(document.sections.iter().any(|s| s.module == "rna"));

        write_report(&document, "report_out").expect("write report");
        assert!(tmp.path().join("report_out.summary.json").exists());
        assert!(tmp.path().join("report_out.report.html").exists());

        env::set_current_dir(cwd).expect("restore cwd");
        fs::remove_dir_all(tmp.path()).ok();
    }

    #[test]
    fn render_html_includes_offline_chart_fallbacks() {
        let document = ReportDocument {
            title: "demo".to_string(),
            generated_by: "test".to_string(),
            sections: vec![ReportSection {
                module: "fastq".to_string(),
                sample: "sample".to_string(),
                source: "/tmp/sample.fastq.summary.json".to_string(),
                metrics: json!({"gc_hist": {"0": 1}, "mean_quality_hist": {"10": 2}}),
            }],
        };

        let html = render_html(&document);
        assert!(html.contains("function drawLineFallback"));
        assert!(html.contains("function drawBarFallback"));
        assert!(html.contains("function drawDoughnutFallback"));
        assert!(html.contains("function renderActiveTabCharts(tabId)"));
    }
}
