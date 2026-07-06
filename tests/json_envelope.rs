use std::io::Write;
use std::process::Command;

use serde_json::Value;

const TABLE: &str = "feature\tS1\tS2\tS3\n\
                     OTU1\t4\t0\t1\n\
                     OTU2\t0\t10\t1\n\
                     OTU3\t3\t10\t1\n\
                     OTU4\t1\t0\t1\n";

fn run_json(table_path: &str) -> Vec<u8> {
    let out = Command::new(env!("CARGO_BIN_EXE_rsomics-beta-diversity"))
        .args([table_path, "-m", "braycurtis", "--json"])
        .output()
        .expect("spawn binary");
    assert!(out.status.success(), "binary exited non-zero: {:?}", out);
    out.stdout
}

#[test]
fn json_is_single_populated_envelope() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("counts.tsv");
    std::fs::File::create(&path)
        .unwrap()
        .write_all(TABLE.as_bytes())
        .unwrap();

    let stdout = run_json(path.to_str().unwrap());

    // A single JSON document — not text + a trailing envelope, and not two docs.
    let doc: Value = serde_json::from_slice(&stdout).unwrap_or_else(|e| {
        panic!(
            "stdout is not one JSON document: {e}\n{}",
            String::from_utf8_lossy(&stdout)
        )
    });

    assert_eq!(doc["status"], "ok");
    let result = &doc["result"];
    assert!(
        !result.is_null(),
        "result must be populated, got null: {doc}"
    );
    assert_eq!(result["metric"], "braycurtis");
    assert_eq!(result["ids"], serde_json::json!(["S1", "S2", "S3"]));
    let rows = result["distances"].as_array().expect("distances array");
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0][0].as_f64().unwrap(), 0.0);
    assert_eq!(rows[1][2], result["distances"][2][1]);
}
