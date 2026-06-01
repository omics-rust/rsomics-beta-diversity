use std::path::PathBuf;
use std::process::Command;

const METRICS: [&str; 5] = [
    "braycurtis",
    "jaccard",
    "euclidean",
    "canberra",
    "cityblock",
];

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-beta-diversity"))
}

fn golden(name: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn oracle_script() -> String {
    format!("{}/tests/oracle_skbio.py", env!("CARGO_MANIFEST_DIR"))
}

/// scikit-bio is the named oracle; skip loudly if it (or python) is unavailable.
fn skbio_python() -> Option<String> {
    for py in ["python3", "python"] {
        let probe = Command::new(py)
            .args(["-c", "import skbio.diversity"])
            .output();
        if let Ok(out) = probe
            && out.status.success()
        {
            return Some(py.to_string());
        }
    }
    eprintln!("SKIP: scikit-bio not importable — install `scikit-bio` to run the differential");
    None
}

fn ours_output(table: &str, metric: &str) -> String {
    let out = Command::new(ours())
        .arg(golden(table))
        .args(["-m", metric])
        .output()
        .expect("run rsomics-beta-diversity");
    assert!(
        out.status.success(),
        "ours failed ({metric}): {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn skbio_output(py: &str, table: &str, metric: &str) -> String {
    let out = Command::new(py)
        .arg(oracle_script())
        .arg(golden(table))
        .arg(metric)
        .output()
        .expect("run scikit-bio oracle");
    assert!(
        out.status.success(),
        "oracle failed ({metric}): {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn differential(table: &str) {
    let Some(py) = skbio_python() else { return };
    for metric in METRICS {
        let ours = ours_output(table, metric);
        let theirs = skbio_output(&py, table, metric);
        assert_eq!(
            ours, theirs,
            "rsomics-beta-diversity output differs from scikit-bio ({table}, {metric})"
        );
    }
}

#[test]
fn matches_skbio_counts_table() {
    differential("counts.tsv");
}

#[test]
fn matches_skbio_edge_table() {
    differential("edge.tsv");
}
