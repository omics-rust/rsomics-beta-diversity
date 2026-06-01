use rsomics_common::{Result, RsomicsError};

/// Between-sample beta-diversity distances, matching the non-phylogenetic
/// metrics `skbio.diversity.beta_diversity` delegates to `scipy.spatial.distance`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Metric {
    BrayCurtis,
    Jaccard,
    Euclidean,
    Canberra,
    Cityblock,
}

impl Metric {
    pub const ALL: [Metric; 5] = [
        Metric::BrayCurtis,
        Metric::Jaccard,
        Metric::Euclidean,
        Metric::Canberra,
        Metric::Cityblock,
    ];

    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Metric::BrayCurtis => "braycurtis",
            Metric::Jaccard => "jaccard",
            Metric::Euclidean => "euclidean",
            Metric::Canberra => "canberra",
            Metric::Cityblock => "cityblock",
        }
    }

    /// # Errors
    /// Errors on an unknown metric token.
    pub fn parse(token: &str) -> Result<Metric> {
        let m = match token {
            "braycurtis" | "bray-curtis" | "bray" => Metric::BrayCurtis,
            "jaccard" => Metric::Jaccard,
            "euclidean" => Metric::Euclidean,
            "canberra" => Metric::Canberra,
            "cityblock" | "manhattan" | "l1" => Metric::Cityblock,
            other => {
                return Err(RsomicsError::InvalidInput(format!(
                    "unknown metric '{other}' (braycurtis|jaccard|euclidean|canberra|cityblock)"
                )));
            }
        };
        Ok(m)
    }

    /// Distance between two equal-length count vectors.
    #[must_use]
    pub fn distance(self, x: &[f64], y: &[f64]) -> f64 {
        match self {
            Metric::BrayCurtis => braycurtis(x, y),
            Metric::Jaccard => jaccard(x, y),
            Metric::Euclidean => euclidean(x, y),
            Metric::Canberra => canberra(x, y),
            Metric::Cityblock => cityblock(x, y),
        }
    }
}

fn braycurtis(x: &[f64], y: &[f64]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        num += (a - b).abs();
        den += a + b;
    }
    num / den
}

/// scipy boolean Jaccard on the count vectors: presence/absence per feature.
/// Distance = mismatched-presence positions / positions present in either.
/// Empty union (both samples empty) yields 0, matching scipy.
fn jaccard(x: &[f64], y: &[f64]) -> f64 {
    let mut neq = 0u64;
    let mut union = 0u64;
    for (&a, &b) in x.iter().zip(y) {
        let pa = a != 0.0;
        let pb = b != 0.0;
        if pa || pb {
            union += 1;
            if pa != pb {
                neq += 1;
            }
        }
    }
    if union == 0 {
        return 0.0;
    }
    neq as f64 / union as f64
}

fn euclidean(x: &[f64], y: &[f64]) -> f64 {
    let mut s = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        let d = a - b;
        s += d * d;
    }
    s.sqrt()
}

/// scipy Canberra: a zero/zero feature contributes nothing. Counts are
/// non-negative, so the denominator is `a + b`; nudging a zero denominator to 1
/// turns the would-be `0/0` into `0/1 = 0` without a per-element branch, which
/// lets the loop autovectorize.
fn canberra(x: &[f64], y: &[f64]) -> f64 {
    let mut s = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        let den = a + b;
        s += (a - b).abs() / (den + f64::from(den == 0.0));
    }
    s
}

fn cityblock(x: &[f64], y: &[f64]) -> f64 {
    let mut s = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        s += (a - b).abs();
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    const A: [f64; 7] = [4.0, 0.0, 3.0, 1.0, 2.0, 0.0, 5.0];
    const B: [f64; 7] = [0.0, 10.0, 10.0, 0.0, 5.0, 1.0, 0.0];

    #[test]
    fn braycurtis_matches_skbio() {
        assert!((Metric::BrayCurtis.distance(&A, &B) - 0.756_097_560_975_609_8).abs() < 1e-15);
    }

    #[test]
    fn jaccard_presence_absence() {
        assert!((Metric::Jaccard.distance(&A, &B) - 0.714_285_714_285_714_3).abs() < 1e-15);
    }

    #[test]
    fn euclidean_matches_skbio() {
        assert!((Metric::Euclidean.distance(&A, &B) - 14.177_446_878_757_825).abs() < 1e-12);
    }

    #[test]
    fn canberra_skips_double_zero() {
        assert!((Metric::Canberra.distance(&A, &B) - 5.967_032_967_032_967).abs() < 1e-12);
    }

    #[test]
    fn cityblock_matches_skbio() {
        assert!((Metric::Cityblock.distance(&A, &B) - 31.0).abs() < 1e-12);
    }

    #[test]
    fn braycurtis_both_empty_is_nan() {
        let z = [0.0; 4];
        assert!(Metric::BrayCurtis.distance(&z, &z).is_nan());
    }

    #[test]
    fn jaccard_empty_union_is_zero() {
        let z = [0.0; 4];
        assert_eq!(Metric::Jaccard.distance(&z, &z), 0.0);
    }
}
