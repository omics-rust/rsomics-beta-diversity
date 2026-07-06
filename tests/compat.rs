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

/// Header-only table (zero features → empty sample vectors). skbio's driver
/// returns an all-zero DistanceMatrix; Bray-Curtis must not leak 0/0 = nan.
#[test]
fn matches_skbio_empty_feature_table() {
    differential("empty.tsv");
}

/// Frozen-oracle goldens for float layouts skbio renders in scientific notation
/// (tiny distances < 1e-4, huge distances >= 1e16). The expected bytes were
/// generated from scikit-bio's own `DistanceMatrix.write`, so this runs offline
/// and pins the Python `repr(float)` formatting without needing skbio present.
fn matches_frozen_golden(table: &str, metric: &str, expected_file: &str) {
    let expected = std::fs::read_to_string(golden(expected_file)).unwrap();
    assert_eq!(
        ours_output(table, metric),
        expected,
        "rsomics-beta-diversity output differs from frozen skbio golden ({expected_file})"
    );
}

#[test]
fn tiny_braycurtis_scientific_notation() {
    matches_frozen_golden("tiny.tsv", "braycurtis", "tiny.braycurtis.expected");
}

#[test]
fn huge_cityblock_scientific_notation() {
    matches_frozen_golden("huge.tsv", "cityblock", "huge.cityblock.expected");
}

#[test]
fn huge_euclidean_scientific_notation() {
    matches_frozen_golden("huge.tsv", "euclidean", "huge.euclidean.expected");
}

#[test]
fn empty_feature_table_is_all_zero() {
    matches_frozen_golden("empty.tsv", "braycurtis", "empty.braycurtis.expected");
}
