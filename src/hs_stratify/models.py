"""Modèles de classification à partir de l'expression des gènes."""

from __future__ import annotations

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


class TopVarianceGenes(BaseEstimator, TransformerMixin):
    """Garde les `k` gènes les plus variables, calculés sur les seules données d'entraînement."""

    def __init__(self, k: int = 2000):
        self.k = k

    def fit(self, X, y=None):
        variances = np.var(np.asarray(X), axis=0)
        self.selected_ = np.sort(np.argsort(variances)[::-1][: self.k])
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.selected_]


def baseline_model(C: float = 0.1, n_genes: int | None = 2000) -> Pipeline:
    """Baseline : régression logistique pénalisée (elastic net), adaptée à p >> n."""
    steps = [] if n_genes is None else [("genes", TopVarianceGenes(n_genes))]
    steps += [
        ("scale", StandardScaler()),
        (
            "clf",
            LogisticRegression(
                solver="saga",
                l1_ratio=0.5,
                C=C,
                class_weight="balanced",
                max_iter=5000,
            ),
        ),
    ]
    return Pipeline(steps)
