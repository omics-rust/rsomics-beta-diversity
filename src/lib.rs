use std::io::{BufRead, Write};

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

mod metric;
pub use metric::Metric;

/// A feature-by-sample count table transposed to per-sample count vectors.
///
/// Input layout (scikit-bio / QIIME / phyloseq convention): the first column is
/// the feature (OTU/taxon) ID, the header row names the samples, cell
/// `[feature][sample]` is the count. Beta-diversity works on the sample vectors,
/// so we store one dense f64 vector per sample.
pub struct CountTable {
    pub sample_names: Vec<String>,
    /// One count vector per sample, each of length `n_features`.
    pub samples: Vec<Vec<f64>>,
}

impl CountTable {
    /// # Errors
    /// Errors on a missing header, a ragged row, a non-numeric cell, or a
    /// negative count (scikit-bio rejects negative counts).
    pub fn parse<R: BufRead>(reader: R, delim: char) -> Result<CountTable> {
        let mut lines = reader.lines();
        let header = loop {
            match lines.next() {
                Some(line) => {
                    let line = line.map_err(RsomicsError::Io)?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    break line;
                }
                None => return Err(RsomicsError::InvalidInput("empty count table".into())),
            }
        };
        let sample_names: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        if sample_names.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "header has no sample columns (need feature-ID column + ≥1 sample)".into(),
            ));
        }
        let n = sample_names.len();
        let mut samples: Vec<Vec<f64>> = vec![Vec::new(); n];
        for (row_idx, line) in lines.enumerate() {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split(delim);
            let feature = fields.next().unwrap_or("");
            let mut seen = 0usize;
            for (col, field) in fields.enumerate() {
                if col >= n {
                    return Err(RsomicsError::InvalidInput(format!(
                        "row {} (feature '{feature}') has more columns than the header",
                        row_idx + 2
                    )));
                }
                let count: f64 = field.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "row {} (feature '{feature}'), sample '{}': '{}' is not a numeric count",
                        row_idx + 2,
                        sample_names[col],
                        field.trim()
                    ))
                })?;
                if count < 0.0 {
                    return Err(RsomicsError::InvalidInput(format!(
                        "row {} (feature '{feature}'), sample '{}': counts cannot be negative",
                        row_idx + 2,
                        sample_names[col]
                    )));
                }
                samples[col].push(count);
                seen += 1;
            }
            if seen != n {
                return Err(RsomicsError::InvalidInput(format!(
                    "row {} (feature '{feature}') has {seen} count columns, header has {n}",
                    row_idx + 2
                )));
            }
        }
        Ok(CountTable {
            sample_names,
            samples,
        })
    }
}

pub struct Config {
    pub metric: Metric,
    pub delim: char,
}

/// A symmetric pairwise distance matrix over the samples, row-major dense.
pub struct DistanceMatrix {
    ids: Vec<String>,
    data: Vec<f64>,
}

impl DistanceMatrix {
    /// Compute the pairwise distance matrix; the upper triangle is evaluated in
    /// parallel over sample pairs, then mirrored.
    #[must_use]
    pub fn compute(table: &CountTable, metric: Metric) -> DistanceMatrix {
        let n = table.samples.len();
        let pairs: Vec<(usize, usize)> = (0..n)
            .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
            .collect();
        let mut data = vec![0.0_f64; n * n];
        let upper: Vec<f64> = pairs
            .par_iter()
            .map(|&(i, j)| metric.distance(&table.samples[i], &table.samples[j]))
            .collect();
        for (&(i, j), &d) in pairs.iter().zip(&upper) {
            data[i * n + j] = d;
            data[j * n + i] = d;
        }
        DistanceMatrix {
            ids: table.sample_names.clone(),
            data,
        }
    }

    /// Write in scikit-bio's `DistanceMatrix` TSV (LSMat) format: an empty
    /// top-left cell then the sample IDs as the header, then one row per sample
    /// (ID + tab-separated distances) with Python `repr(float)` formatting.
    ///
    /// # Errors
    /// Propagates write errors.
    pub fn write_tsv<W: Write>(&self, mut out: W) -> Result<()> {
        let n = self.ids.len();
        for id in &self.ids {
            write!(out, "\t{id}").map_err(RsomicsError::Io)?;
        }
        writeln!(out).map_err(RsomicsError::Io)?;
        let mut row = String::new();
        for i in 0..n {
            row.clear();
            row.push_str(&self.ids[i]);
            for j in 0..n {
                row.push('\t');
                push_pyrepr(&mut row, self.data[i * n + j]);
            }
            writeln!(out, "{row}").map_err(RsomicsError::Io)?;
        }
        Ok(())
    }
}

/// Append `x` formatted as Python's `repr(float)` would: shortest round-trip
/// decimal, but integer-valued finite floats keep a trailing `.0`, and NaN
/// renders lowercase `nan` (skbio emits `nan` for Bray-Curtis of two empty
/// samples).
fn push_pyrepr(buf: &mut String, x: f64) {
    use std::fmt::Write;
    if x.is_nan() {
        buf.push_str("nan");
        return;
    }
    if x.is_infinite() {
        buf.push_str(if x < 0.0 { "-inf" } else { "inf" });
        return;
    }
    let start = buf.len();
    let _ = write!(buf, "{x}");
    if !buf[start..].contains(['.', 'e', 'E']) {
        buf.push_str(".0");
    }
}

/// # Errors
/// Propagates parse and write errors.
pub fn run<R: BufRead, W: Write>(reader: R, out: W, cfg: &Config) -> Result<()> {
    let table = CountTable::parse(reader, cfg.delim)?;
    let dm = DistanceMatrix::compute(&table, cfg.metric);
    dm.write_tsv(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table() -> &'static str {
        "feature\tS1\tS2\tS3\n\
         OTU1\t4\t0\t1\n\
         OTU2\t0\t10\t1\n\
         OTU3\t3\t10\t1\n\
         OTU4\t1\t0\t1\n\
         OTU5\t2\t5\t1\n\
         OTU6\t0\t1\t1\n\
         OTU7\t5\t0\t1\n"
    }

    #[test]
    fn parses_transposed_columns() {
        let t = CountTable::parse(table().as_bytes(), '\t').unwrap();
        assert_eq!(t.sample_names, ["S1", "S2", "S3"]);
        assert_eq!(t.samples[0], [4.0, 0.0, 3.0, 1.0, 2.0, 0.0, 5.0]);
    }

    #[test]
    fn braycurtis_matrix_matches_skbio() {
        let t = CountTable::parse(table().as_bytes(), '\t').unwrap();
        let dm = DistanceMatrix::compute(&t, Metric::BrayCurtis);
        let mut buf = Vec::new();
        dm.write_tsv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert_eq!(
            s,
            "\tS1\tS2\tS3\n\
             S1\t0.0\t0.7560975609756098\t0.5454545454545454\n\
             S2\t0.7560975609756098\t0.0\t0.7575757575757576\n\
             S3\t0.5454545454545454\t0.7575757575757576\t0.0\n"
        );
    }

    #[test]
    fn cityblock_integers_keep_dot_zero() {
        let t = CountTable::parse(table().as_bytes(), '\t').unwrap();
        let dm = DistanceMatrix::compute(&t, Metric::Cityblock);
        let mut buf = Vec::new();
        dm.write_tsv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("S1\t0.0\t31.0\t12.0\n"));
    }

    #[test]
    fn negative_count_errors() {
        let bad = "feature\tA\tB\nOTU1\t1\t-2\n";
        assert!(CountTable::parse(bad.as_bytes(), '\t').is_err());
    }

    #[test]
    fn ragged_row_errors() {
        let bad = "feature\tA\tB\nOTU1\t4\n";
        assert!(CountTable::parse(bad.as_bytes(), '\t').is_err());
    }
}
