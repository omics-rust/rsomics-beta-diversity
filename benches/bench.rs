use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

use criterion::{Criterion, criterion_group, criterion_main};

fn bench_beta(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-beta-diversity");
    let table = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/counts.tsv");
    c.bench_function("rsomics-beta-diversity braycurtis golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .arg(&table)
                .args(["-m", "braycurtis", "-t", "1"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_beta);
criterion_main!(benches);
