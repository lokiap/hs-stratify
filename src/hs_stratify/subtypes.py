"""Sous-types moléculaires des lésions : clustering de consensus et stabilité entre cohortes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans


@dataclass
class ConsensusResult:
    k: int
    labels: np.ndarray
    consensus: np.ndarray
    pac: float


def consensus_matrix(
    X: np.ndarray, k: int, n_resamples: int = 100, fraction: float = 0.8, seed: int = 0
) -> np.ndarray:
    """Fréquence à laquelle deux échantillons tombent dans le même groupe (Monti et al., 2003).

    À chaque tirage, on garde `fraction` des échantillons et on lance un k-means. Le consensus
    est le nombre de fois où deux échantillons sont regroupés, divisé par le nombre de fois où
    ils ont été tirés ensemble.
    """
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    together = np.zeros((n, n))
    sampled = np.zeros((n, n))
    for i in range(n_resamples):
        idx = np.sort(rng.choice(n, size=int(round(fraction * n)), replace=False))
        labels = KMeans(n_clusters=k, n_init=5, random_state=seed + i).fit_predict(X[idx])
        same = (labels[:, None] == labels[None, :]).astype(float)
        together[np.ix_(idx, idx)] += same
        sampled[np.ix_(idx, idx)] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        consensus = np.where(sampled > 0, together / sampled, 0.0)
    np.fill_diagonal(consensus, 1.0)
    return consensus


def pac_score(consensus: np.ndarray, low: float = 0.1, high: float = 0.9) -> float:
    """Part des paires ambiguës (consensus entre `low` et `high`) : plus bas = plus net."""
    values = consensus[np.triu_indices_from(consensus, k=1)]
    return float(np.mean((values > low) & (values < high)))


def consensus_clustering(
    X: np.ndarray, k_values: tuple[int, ...] = (2, 3, 4, 5), n_resamples: int = 100, seed: int = 0
) -> dict[int, ConsensusResult]:
    """Consensus pour chaque k, puis groupes finaux par classification hiérarchique du consensus."""
    results = {}
    for k in k_values:
        consensus = consensus_matrix(X, k, n_resamples=n_resamples, seed=seed)
        distance = squareform(1.0 - consensus, checks=False)
        tree = linkage(distance, method="average")
        labels = fcluster(tree, t=k, criterion="maxclust") - 1
        results[k] = ConsensusResult(k, labels, consensus, pac_score(consensus))
    return results


def best_k(results: dict[int, ConsensusResult]) -> int:
    """Le k au PAC le plus bas ; à égalité, le plus petit k."""
    return min(results, key=lambda k: (round(results[k].pac, 3), k))


def nearest_centroid(X_train: np.ndarray, labels: np.ndarray, X_new: np.ndarray) -> np.ndarray:
    """Attribue chaque nouvel échantillon au sous-type dont le profil moyen est le plus proche."""
    classes = np.unique(labels)
    centroids = np.vstack([X_train[labels == c].mean(axis=0) for c in classes])
    distances = ((X_new[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    return classes[np.argmin(distances, axis=1)]


def align_labels(labels: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Renumérote `labels` pour qu'ils correspondent au mieux à `reference` (même k).

    Les numéros de groupe d'un clustering sont arbitraires : sans cet alignement, deux
    partitions identiques peuvent porter des numéros échangés.
    """
    classes = np.unique(np.r_[labels, reference])
    overlap = np.array(
        [[np.sum((labels == a) & (reference == b)) for b in classes] for a in classes]
    )
    rows, cols = linear_sum_assignment(-overlap)
    mapping = {classes[r]: classes[c] for r, c in zip(rows, cols, strict=True)}
    return np.array([mapping[v] for v in labels])


def marker_genes(X: pd.DataFrame, labels: np.ndarray, n: int = 15) -> pd.DataFrame:
    """Gènes les plus surexprimés dans chaque sous-type par rapport aux autres."""
    rows = []
    for c in np.unique(labels):
        inside = X[labels == c].mean()
        outside = X[labels != c].mean()
        diff = (inside - outside).sort_values(ascending=False).head(n)
        rows += [{"subtype": int(c), "gene": g, "difference": round(v, 3)} for g, v in diff.items()]
    return pd.DataFrame(rows)
