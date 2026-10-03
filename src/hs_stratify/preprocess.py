"""Harmonisation des matrices d'expression entre cohortes."""

from __future__ import annotations

import pandas as pd


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
