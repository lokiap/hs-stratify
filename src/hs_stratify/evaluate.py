"""Évaluation croisée par cohorte."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd
from sklearn.base import ClassifierMixin
from sklearn.metrics import roc_auc_score

from hs_stratify.splits import leave_one_cohort_out


def evaluate_loco(
    X: np.ndarray,
    y: np.ndarray,
    cohorts: Sequence[str],
    make_model: Callable[[], ClassifierMixin],
) -> pd.DataFrame:
    """Entraîne sur toutes les cohortes sauf une et mesure l'AUROC sur celle laissée de côté."""
    rows = []
    for held_out, train, test in leave_one_cohort_out(cohorts):
        y_test = y[test]
        if len(np.unique(y_test)) < 2:
            rows.append({"cohort": held_out, "n_test": len(test), "auroc": np.nan})
            continue
        model = make_model().fit(X[train], y[train])
        scores = model.predict_proba(X[test])[:, 1]
        rows.append(
            {"cohort": held_out, "n_test": len(test), "auroc": roc_auc_score(y_test, scores)}
        )
    return pd.DataFrame(rows)
