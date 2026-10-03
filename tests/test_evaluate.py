import numpy as np

from hs_stratify.evaluate import evaluate_loco
from hs_stratify.models import baseline_model


def test_loco_on_synthetic_signal():
    rng = np.random.default_rng(42)
    n, p = 90, 50
    y = rng.integers(0, 2, size=n)
    X = rng.normal(size=(n, p))
    X[:, 0] += 2.0 * y  # un gène porte le signal
    cohorts = np.repeat(["A", "B", "C"], n // 3)

    results = evaluate_loco(X, y, cohorts, lambda: baseline_model(C=1.0))

    assert list(results["cohort"]) == ["A", "B", "C"]
    assert (results["auroc"] > 0.7).all()
