import numpy as np
import pandas as pd

from hs_stratify.preprocess import (
    common_genes,
    ensure_log2,
    log_cpm,
    probes_to_genes,
    standardize_per_cohort,
)


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


def test_log_cpm_filters_low_genes():
    counts = pd.DataFrame(
        {"S1": [1_000_000, 0], "S2": [2_000_000, 1]}, index=["KRT16", "RARE"], dtype=float
    )
    result = log_cpm(counts, min_cpm=1.0, min_fraction=0.6)
    assert list(result.index) == ["KRT16"]
    np.testing.assert_allclose(result.loc["KRT16"], np.log2(1e6 + 1))


def test_ensure_log2_only_transforms_raw_intensities():
    already_log = pd.DataFrame({"S1": [7.0, 12.0]})
    assert ensure_log2(already_log).equals(already_log)
    raw = pd.DataFrame({"S1": [1023.0, 3.0]})
    np.testing.assert_allclose(ensure_log2(raw)["S1"], [10.0, 2.0])


def test_probes_to_genes_averages_probes():
    values = pd.DataFrame({"S1": [2.0, 4.0, 9.0]}, index=["p1", "p2", "p3"])
    mapping = pd.Series({"p1": "TNF", "p2": "TNF"})
    result = probes_to_genes(values, mapping)
    assert result.to_dict() == {"S1": {"TNF": 3.0}}
