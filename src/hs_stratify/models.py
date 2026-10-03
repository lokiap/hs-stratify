"""Modèles de prédiction de la réponse au traitement."""

from __future__ import annotations

from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def baseline_model(C: float = 0.1) -> Pipeline:
    """Baseline : régression logistique pénalisée (elastic net), adaptée à p >> n."""
    return Pipeline(
        [
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
    )
