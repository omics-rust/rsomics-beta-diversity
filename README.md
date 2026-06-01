# rsomics-beta-diversity

Pairwise **between-sample** beta-diversity distance matrix from a feature
(OTU/ASV/taxon) count table.

Reads a feature-by-sample TSV — first column is the feature ID, the header names
the samples, cell `[feature][sample]` is a count — and writes a square symmetric
distance matrix in scikit-bio's `DistanceMatrix` TSV format.

```
rsomics-beta-diversity counts.tsv
cat counts.tsv | rsomics-beta-diversity -m jaccard -t 8
```

## Metrics

| `-m` value            | distance |
|-----------------------|----------|
| `braycurtis`          | `Σ\|xᵢ-yᵢ\| / Σ(xᵢ+yᵢ)` — the microbiome default (`nan` if both samples are empty) |
| `jaccard`             | presence/absence: mismatched features / features present in either sample |
| `euclidean`           | `√Σ(xᵢ-yᵢ)²` |
| `canberra`            | `Σ \|xᵢ-yᵢ\| / (\|xᵢ\|+\|yᵢ\|)`, dropping features absent from both |
| `cityblock` (`manhattan`) | `Σ\|xᵢ-yᵢ\|` |

These are exactly the non-phylogenetic metrics
`skbio.diversity.beta_diversity` delegates to `scipy.spatial.distance`. Negative
counts are rejected, matching scikit-bio.

### Output format

scikit-bio's `DistanceMatrix` TSV: an empty top-left cell then the sample IDs as
the header row, then one row per sample (ID followed by tab-separated
distances). Floats use Python's shortest round-trip `repr`, so the output is
byte-identical to `skbio.DistanceMatrix.write()`.

### Not included: UniFrac

Unweighted and weighted **UniFrac** are phylogenetic beta-diversity metrics —
they need a rooted phylogenetic tree relating the features, not just the count
table. That is a different operation with a different input contract and lives
in its own crate; it is intentionally out of scope here.

## Origin

This crate is an independent Rust reimplementation of the non-phylogenetic
beta-diversity operation provided by `scikit-bio`
(`skbio.diversity.beta_diversity`, which delegates to
`scipy.spatial.distance.pdist`), based on:

- The published distance definitions (Bray & Curtis 1957; Jaccard 1912; the
  Canberra metric of Lance & Williams 1966; Euclidean and Manhattan/cityblock).
- The black-box behaviour of `skbio.diversity.beta_diversity`: presence/absence
  Jaccard on the count vectors, the double-zero exclusion in Canberra, the
  empty-sample conventions (`braycurtis` → `nan`, `jaccard` → `0`), and the
  `DistanceMatrix` TSV serialization including float `repr`.

scikit-bio is BSD-3-Clause and was read and cited. Test fixtures are
independently generated count tables.

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio <https://scikit-bio.org> (BSD-3-Clause).

## Compatibility & performance

`tests/compat.rs` runs this binary and the scikit-bio oracle
(`tests/oracle_skbio.py`) for every metric on the golden tables and asserts
byte-identical output. The differential is skipped loudly when scikit-bio is not
importable.

The O(samples²·features) pair loop is parallelised with rayon over the matrix's
upper triangle. On a 2000-sample × 5000-feature table, single-threaded vs
scikit-bio 0.7.2 (scipy `pdist`): jaccard 3.0×, canberra 1.28×, braycurtis
1.20×, euclidean/cityblock ~1.0× (near parity — those two scipy primitives are
the most tightly tuned). The pair loop adds ~4.9× more on 8 threads, which
scikit-bio cannot do for these metrics.
