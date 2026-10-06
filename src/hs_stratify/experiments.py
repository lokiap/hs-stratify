"""Expériences reproductibles à partir des matrices construites par `hs-stratify build`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from hs_stratify.cohorts import Cohort
from hs_stratify.evaluate import evaluate_loco
from hs_stratify.labels import tissue_labels
from hs_stratify.models import baseline_model
from hs_stratify.preprocess import common_genes, standardize_per_cohort, standardize_per_sample

# Lésion contre peau sans lésion (non lésionnelle chez un patient, ou peau saine).
# La peau péri-lésionnelle, entre les deux, est écartée.
LESIONAL_TARGET = {"lesional": 1, "nonlesional": 0, "healthy": 0}


def load_processed(cohort_id: str, processed_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.read_csv(processed_dir / f"{cohort_id}_expression.csv.gz", index_col=0)
    meta = pd.read_csv(processed_dir / f"{cohort_id}_metadata.csv.gz", index_col=0, dtype=str)
    return expression, meta


def lesional_dataset(
    data: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    cohorts: list[Cohort],
    normalization: str = "sample",
) -> tuple[pd.DataFrame, pd.Series]:
    """Assemble X (échantillons x gènes communs) et y (1 = lésion), indexés par cohorte."""
    matrices, targets = {}, {}
    for cohort in cohorts:
        if cohort.tissue is None or cohort.id not in data:
            continue
        expression, meta = data[cohort.id]
        tissue = tissue_labels(meta, cohort.tissue)
        y = tissue.map(LESIONAL_TARGET).dropna().astype(int)
        y = y[y.index.isin(expression.index)]
        if y.empty:
            continue
        matrices[cohort.id] = expression.loc[y.index]
        targets[cohort.id] = y

    genes = common_genes(matrices)
    if normalization == "cohort":
        X = standardize_per_cohort(matrices)
    elif normalization == "sample":
        X = pd.concat(
            {cid: standardize_per_sample(m[genes]) for cid, m in matrices.items()},
            names=["cohort", "sample"],
        )
    else:
        raise ValueError(f"Normalisation inconnue : {normalization}")
    y = pd.concat(targets, names=["cohort", "sample"]).loc[X.index]
    return X, y


def run_lesional(
    data: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    cohorts: list[Cohort],
    normalizations: tuple[str, ...] = ("sample", "cohort"),
) -> pd.DataFrame:
    """AUROC lésion / non-lésion sur chaque cohorte laissée de côté, pour chaque normalisation."""
    tables = []
    for normalization in normalizations:
        X, y = lesional_dataset(data, cohorts, normalization)
        groups = X.index.get_level_values("cohort")
        table = evaluate_loco(X.to_numpy(), y.to_numpy(), groups, baseline_model)
        counts = y.groupby(level="cohort").agg(n_lesional="sum", n_total="size")
        table = table.join(counts, on="cohort")
        table.insert(0, "normalization", normalization)
        table["n_genes"] = X.shape[1]
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def summarize(results: pd.DataFrame) -> pd.DataFrame:
    """Une ligne par normalisation : AUROC sur les cohortes qui ont les deux classes."""
    valid = results.dropna(subset=["auroc"])
    return valid.groupby("normalization")["auroc"].agg(["size", "mean", "median", "min"]).round(3)
