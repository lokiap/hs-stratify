import numpy as np
import pandas as pd

from hs_stratify.cohorts import Cohort
from hs_stratify.experiments import lesional_dataset, run_lesional, summarize

TISSUE = {"field": "ch1:tissue", "classes": {"lesional": "LS", "nonlesional": "NL"}}


def _cohort(rng, n, offset, only_lesional=False):
    genes = [f"G{i}" for i in range(30)]
    tissue = ["LS"] * n if only_lesional else ["LS", "NL"] * (n // 2)
    is_lesion = np.array([t == "LS" for t in tissue])
    values = rng.normal(5.0 + offset, 1.0, size=(len(tissue), len(genes)))
    values[:, 0] += 3.0 * is_lesion  # un gène marqueur de lésion
    values[:, 1] -= 3.0 * is_lesion
    index = [f"S{offset}_{i}" for i in range(len(tissue))]
    expression = pd.DataFrame(values, index=index, columns=genes)
    meta = pd.DataFrame({"ch1:tissue": tissue}, index=index)
    return expression, meta


def test_lesional_dataset_shapes_and_labels():
    rng = np.random.default_rng(0)
    data = {"A": _cohort(rng, 20, 0), "B": _cohort(rng, 10, 4, only_lesional=True)}
    cohorts = [Cohort(cid, "discovery", "rnaseq", tissue=TISSUE) for cid in data]

    X, y = lesional_dataset(data, cohorts, "sample")

    assert X.shape == (30, 30)
    assert list(X.index.names) == ["cohort", "sample"]
    assert y.loc["B"].eq(1).all()
    assert y.loc["A"].sum() == 10


def test_run_lesional_generalizes_across_cohorts():
    rng = np.random.default_rng(1)
    data = {cid: _cohort(rng, 30, off) for cid, off in [("A", 0), ("B", 3), ("C", -2)]}
    cohorts = [Cohort(cid, "discovery", "rnaseq", tissue=TISSUE) for cid in data]

    results = run_lesional(data, cohorts)

    assert set(results["normalization"]) == {"sample", "cohort"}
    assert (results["auroc"] > 0.9).all()
    assert set(summarize(results).index) == {"sample", "cohort"}
