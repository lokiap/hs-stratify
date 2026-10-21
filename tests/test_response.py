import numpy as np
import pandas as pd

from hs_stratify.experiments import run_response
from hs_stratify.models import ridge_model as baseline_model
from hs_stratify.response import (
    bootstrap_auroc_ci,
    leave_one_out_scores,
    permutation_p_value,
)


def make_model():
    return baseline_model(C=1.0, n_genes=None)


def _patients(rng, n, signal):
    y = np.repeat([1, 0], n // 2)
    X = rng.normal(size=(n, 20))
    X[:, 0] += signal * y
    return X, y


def test_bootstrap_ci_contains_observed_auroc():
    y = np.array([1, 1, 1, 0, 0, 0])
    scores = np.array([0.9, 0.8, 0.3, 0.4, 0.2, 0.1])
    low, high = bootstrap_auroc_ci(y, scores, n_bootstrap=200)
    assert low <= 8 / 9 <= high


def test_permutation_separates_signal_from_noise():
    rng = np.random.default_rng(0)
    X_tr, y_tr = _patients(rng, 20, signal=4.0)
    X_te, y_te = _patients(rng, 12, signal=4.0)
    p_signal = permutation_p_value(make_model, X_tr, y_tr, X_te, y_te, 1.0, n_permutations=30)
    assert p_signal < 0.1

    X_tr, y_tr = _patients(rng, 20, signal=0.0)
    X_te, y_te = _patients(rng, 12, signal=0.0)
    p_noise = permutation_p_value(make_model, X_tr, y_tr, X_te, y_te, 0.5, n_permutations=30)
    assert p_noise > 0.1


def test_leave_one_out_never_scores_a_patient_it_trained_on():
    rng = np.random.default_rng(1)
    X, y = _patients(rng, 10, signal=3.0)
    scores = leave_one_out_scores(make_model, X, y)
    assert scores.shape == (10,)
    assert ((scores >= 0) & (scores <= 1)).all()


def test_run_response_reports_both_directions_and_pooled():
    rng = np.random.default_rng(2)
    genes = [f"G{i}" for i in range(25)]
    data, responses = {}, {}
    for cid, n in [("A", 12), ("B", 14)]:
        X, y = _patients(rng, n, signal=6.0)
        X = np.c_[X, rng.normal(size=(n, 5))] + 5.0
        index = [f"{cid}{i}" for i in range(n)]
        data[cid] = (pd.DataFrame(X, index=index, columns=genes), pd.DataFrame(index=index))
        responses[cid] = pd.Series(y, index=index)

    table = run_response(
        data,
        responses,
        n_genes=None,
        C=1.0,
        n_permutations=10,
        n_bootstrap=50,
        n_loo_permutations=5,
    )

    assert table["setting"].tolist() == ["cross-cohort", "cross-cohort", "pooled leave-one-out"]
    assert table.loc[0, "n_test"] == "7R/7NR"
    assert (table["auroc"] > 0.7).all()
