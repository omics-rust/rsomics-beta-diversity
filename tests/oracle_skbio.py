#!/usr/bin/env python3
"""scikit-bio beta-diversity oracle for rsomics-beta-diversity compat tests.

Reads a feature-by-sample TSV count table (feature-ID column + sample columns)
on argv[1], runs `skbio.diversity.beta_diversity(metric, ...)` for the metric on
argv[2], and writes the resulting DistanceMatrix in scikit-bio's own TSV form —
the exact bytes rsomics-beta-diversity targets.
"""

import sys

import numpy as np
from skbio.diversity import beta_diversity


def main():
    path = sys.argv[1]
    metric = sys.argv[2]
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith("#")]
    header = lines[0].split("\t")
    ids = header[1:]
    cols = [[] for _ in ids]
    for ln in lines[1:]:
        parts = ln.split("\t")
        for i, v in enumerate(parts[1:]):
            cols[i].append(float(v))
    counts = np.array(cols, dtype=float)
    dm = beta_diversity(metric, counts, ids)
    dm.write(sys.stdout)


if __name__ == "__main__":
    main()
