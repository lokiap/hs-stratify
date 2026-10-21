"""Évaluation de la prédiction de réponse sur de très petits effectifs.

Avec une dizaine de patients par cohorte, une AUROC seule ne veut rien dire : on lui
adjoint un intervalle de confiance (bootstrap) et un test de permutation, qui dit à quelle
fréquence un modèle entraîné sur des étiquettes mélangées fait aussi bien.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from sklearn.base import ClassifierMixin
from sklearn.metrics import roc_auc_score

ModelFactory = Callable[[], ClassifierMixin]


def fit_score(make_model: ModelFactory, X_train, y_train, X_test) -> np.ndarray:
    """Entraîne un modèle neuf et renvoie la probabilité de réponse sur `X_test`."""
    return make_model().fit(X_train, y_train).predict_proba(X_test)[:, 1]


def bootstrap_auroc_ci(
    y: np.ndarray, scores: np.ndarray, n_bootstrap: int = 1000, seed: int = 0
) -> tuple[float, float]:
    """Intervalle à 95 % de l'AUROC, en rééchantillonnant répondeurs et non-répondeurs à part."""
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    values = []
    for _ in range(n_bootstrap):
        idx = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        values.append(roc_auc_score(y[idx], scores[idx]))
    low, high = np.percentile(values, [2.5, 97.5])
    return float(low), float(high)


def permutation_p_value(
    make_model: ModelFactory,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    observed: float,
    n_permutations: int = 200,
    seed: int = 0,
) -> float:
    """Part des modèles entraînés sur des étiquettes mélangées qui atteignent l'AUROC observée."""
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_permutations):
        shuffled = rng.permutation(y_train)
        if roc_auc_score(y_test, fit_score(make_model, X_train, shuffled, X_test)) >= observed:
            hits += 1
    return (hits + 1) / (n_permutations + 1)


def leave_one_out_scores(make_model: ModelFactory, X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Score de chaque patient par un modèle entraîné sur tous les autres."""
    scores = np.empty(len(y))
    for i in range(len(y)):
        train = np.arange(len(y)) != i
        scores[i] = fit_score(make_model, X[train], y[train], X[i : i + 1])[0]
    return scores


def leave_one_out_p_value(
    make_model: ModelFactory,
    X: np.ndarray,
    y: np.ndarray,
    observed: float,
    n_permutations: int = 100,
    seed: int = 0,
) -> float:
    """Test de permutation pour l'AUROC en leave-one-out (un mélange des étiquettes par tour)."""
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_permutations):
        shuffled = rng.permutation(y)
        if roc_auc_score(shuffled, leave_one_out_scores(make_model, X, shuffled)) >= observed:
            hits += 1
    return (hits + 1) / (n_permutations + 1)
