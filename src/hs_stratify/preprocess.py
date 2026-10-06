"""Harmonisation des matrices d'expression entre cohortes."""

from __future__ import annotations

import numpy as np
import pandas as pd


def log_cpm(counts: pd.DataFrame, min_cpm: float = 1.0, min_fraction: float = 0.2) -> pd.DataFrame:
    """Comptages bruts (gènes x échantillons) -> log2(CPM + 1).

    Garde les gènes avec un CPM >= `min_cpm` dans au moins `min_fraction` des échantillons,
    le même filtre que l'étude de référence.
    """
    cpm = counts / counts.sum(axis=0) * 1e6
    keep = (cpm >= min_cpm).mean(axis=1) >= min_fraction
    return np.log2(cpm[keep] + 1)


def ensure_log2(values: pd.DataFrame) -> pd.DataFrame:
    """Passe des intensités de puce en log2 si elles ne le sont pas déjà."""
    if values.max().max() > 100:
        return np.log2(values.clip(lower=0) + 1)
    return values


def probes_to_genes(values: pd.DataFrame, probe_to_symbol: pd.Series) -> pd.DataFrame:
    """Sondes x échantillons -> gènes x échantillons, en moyennant les sondes d'un même gène."""
    mapped = values.join(probe_to_symbol.rename("symbol"), how="inner")
    return mapped.groupby("symbol").mean()


def common_genes(matrices: dict[str, pd.DataFrame]) -> list[str]:
    """Gènes présents dans toutes les cohortes (matrices échantillons x gènes)."""
    genes: set[str] | None = None
    for matrix in matrices.values():
        genes = set(matrix.columns) if genes is None else genes & set(matrix.columns)
    return sorted(genes or [])


def standardize_per_cohort(matrices: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Restreint aux gènes communs et centre-réduit chaque gène au sein de sa cohorte.

    Correction simple de l'effet d'étude : chaque cohorte est ramenée à moyenne 0
    et écart-type 1 par gène, sans utiliser d'information d'une autre cohorte.
    Renvoie une matrice concaténée avec un index (cohorte, échantillon).
    """
    genes = common_genes(matrices)
    if not genes:
        raise ValueError("Aucun gène commun entre les cohortes")

    blocks = []
    for cohort_id, matrix in matrices.items():
        block = matrix[genes]
        std = block.std(ddof=0).replace(0, 1.0)
        block = (block - block.mean()) / std
        block.index = pd.MultiIndex.from_product(
            [[cohort_id], block.index], names=["cohort", "sample"]
        )
        blocks.append(block)
    return pd.concat(blocks)


def standardize_per_sample(matrix: pd.DataFrame) -> pd.DataFrame:
    """Centre-réduit chaque échantillon sur l'ensemble de ses gènes.

    Contrairement à la standardisation par cohorte, le résultat ne dépend pas des autres
    échantillons : une cohorte qui ne contient que des lésions garde son signal.
    """
    std = matrix.std(axis=1, ddof=0).replace(0, 1.0)
    return matrix.sub(matrix.mean(axis=1), axis=0).div(std, axis=0)
