import numpy as np
import pandas as pd

from hs_stratify.preprocess import common_genes, standardize_per_cohort


def _matrix(genes, n, offset):
    rng = np.random.default_rng(0)
    return pd.DataFrame(rng.normal(offset, 2.0, size=(n, len(genes))), columns=genes)


def test_common_genes_intersection():
    matrices = {"A": _matrix(["IL17A", "TNF", "KRT16"], 3, 0), "B": _matrix(["TNF", "KRT16"], 3, 0)}
    assert common_genes(matrices) == ["KRT16", "TNF"]


def test_standardize_removes_cohort_offset():
    matrices = {"A": _matrix(["TNF", "KRT16"], 20, 0.0), "B": _matrix(["TNF", "KRT16"], 20, 50.0)}
    merged = standardize_per_cohort(matrices)

    assert merged.shape == (40, 2)
    for cohort in ["A", "B"]:
        block = merged.xs(cohort, level="cohort")
        np.testing.assert_allclose(block.mean(), 0.0, atol=1e-9)
        np.testing.assert_allclose(block.std(ddof=0), 1.0, atol=1e-9)
