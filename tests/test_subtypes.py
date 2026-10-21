import numpy as np
import pandas as pd

from hs_stratify.cohorts import Cohort
from hs_stratify.experiments import lesional_matrix, run_subtypes
from hs_stratify.subtypes import (
    best_k,
    consensus_clustering,
    nearest_centroid,
    pac_score,
)

TISSUE = {"field": "ch1:tissue", "classes": {"lesional": "LS", "nonlesional": "NL"}}


def _two_groups(rng, n=30, p=40):
    truth = np.repeat([0, 1], n // 2)
    X = rng.normal(size=(n, p))
    X[truth == 1, :10] += 4.0
    return X, truth


def test_pac_is_zero_for_a_perfect_consensus():
    consensus = np.kron(np.eye(2), np.ones((3, 3)))
    assert pac_score(consensus) == 0.0


def test_consensus_finds_two_clear_groups():
    X, truth = _two_groups(np.random.default_rng(0))
    results = consensus_clustering(X, (2, 3, 4), n_resamples=20)
    assert best_k(results) == 2
    labels = results[2].labels
    assert len(set(zip(labels, truth, strict=True))) == 2  # même partition, à l'étiquette près


def test_nearest_centroid_assigns_new_samples():
    X, truth = _two_groups(np.random.default_rng(1))
    X_new, truth_new = _two_groups(np.random.default_rng(2), n=10)
    predicted = nearest_centroid(X, truth, X_new)
    assert (predicted == truth_new).all()


def _cohort(rng, cid, offset, n=24):
    genes = [f"G{i}" for i in range(60)]
    truth = np.repeat([0, 1], n // 2)
    values = rng.normal(5.0 + offset, 1.0, size=(n + 4, len(genes)))
    values[:n][truth == 1, :15] += 4.0
    index = [f"{cid}_{i}" for i in range(n + 4)]
    tissue = ["LS"] * n + ["NL"] * 4
    expression = pd.DataFrame(values, index=index, columns=genes)
    meta = pd.DataFrame({"ch1:tissue": tissue}, index=index)
    response = pd.Series(np.r_[np.ones(n // 2), np.zeros(n // 2)][:n].astype(int), index=index[:n])
    return expression, meta, response


def test_lesional_matrix_keeps_lesions_and_centers_each_cohort():
    rng = np.random.default_rng(3)
    data = {cid: _cohort(rng, cid, off)[:2] for cid, off in [("A", 0), ("B", 10)]}
    cohorts = [Cohort(cid, "discovery", "rnaseq", tissue=TISSUE) for cid in data]
    lesions = lesional_matrix(data, cohorts)
    assert lesions.shape == (48, 60)
    np.testing.assert_allclose(lesions.loc["B"].mean(), 0.0, atol=1e-9)


def test_run_subtypes_recovers_groups_across_cohorts():
    rng = np.random.default_rng(4)
    raw = {cid: _cohort(rng, cid, off) for cid, off in [("A", 0), ("B", 3), ("C", -2)]}
    data = {cid: (e, m) for cid, (e, m, _) in raw.items()}
    responses = {"A": raw["A"][2]}
    cohorts = [Cohort(cid, "discovery", "rnaseq", tissue=TISSUE) for cid in data]

    tables = run_subtypes(data, cohorts, responses, k_values=(2, 3), n_genes=30, n_resamples=20)

    assert tables["pac"].set_index("k")["pac"].idxmin() == 2
    assert (tables["stability"]["ari"] > 0.9).all()
    assert (tables["stability"]["agreement"] > 0.95).all()
    summary = tables["response_summary"]
    assert sorted(summary["rate"]) == [0.0, 1.0]
    assert "fisher_p" in summary.columns
    assert set(tables["markers"]["subtype"]) == {0, 1}


def test_align_labels_swaps_permuted_ids():
    from hs_stratify.subtypes import align_labels

    reference = np.array([0, 0, 1, 1, 2])
    assert align_labels(np.array([2, 2, 0, 0, 1]), reference).tolist() == reference.tolist()
