"""Expériences reproductibles à partir des matrices construites par `hs-stratify build`."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd

from hs_stratify.cohorts import Cohort
from hs_stratify.evaluate import evaluate_loco
from hs_stratify.labels import tissue_labels
from hs_stratify.models import baseline_model, ridge_model
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
    log: Callable[[str], None] | None = None,
) -> pd.DataFrame:
    """AUROC lésion / non-lésion sur chaque cohorte laissée de côté, pour chaque normalisation."""
    tables = []
    for normalization in normalizations:
        X, y = lesional_dataset(data, cohorts, normalization)
        if log:
            log(f"Normalisation {normalization} : {X.shape[0]} échantillons x {X.shape[1]} gènes")
        groups = X.index.get_level_values("cohort")
        table = evaluate_loco(X.to_numpy(), y.to_numpy(), groups, baseline_model, log=log)
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


def load_response(cohort_id: str, processed_dir: Path) -> pd.Series | None:
    """Étiquettes de réponse écrites par `hs-stratify build`, si la cohorte en a."""
    path = processed_dir / f"{cohort_id}_labels.csv"
    if not path.exists():
        return None
    return pd.read_csv(path, index_col=0)["response"]


def lesional_matrix(
    data: dict[str, tuple[pd.DataFrame, pd.DataFrame]], cohorts: list[Cohort]
) -> pd.DataFrame:
    """Biopsies de lésion avant traitement, gènes communs, centrées au sein de chaque cohorte.

    Le centrage retire le niveau propre à chaque étude (technologie, protocole) mais garde les
    écarts entre lésions d'une même étude, qui sont ce que les sous-types doivent capturer.
    """
    blocks = {}
    for cohort in cohorts:
        if cohort.tissue is None or cohort.id not in data:
            continue
        expression, meta = data[cohort.id]
        tissue = tissue_labels(meta, cohort.tissue)
        samples = tissue.index[(tissue == "lesional") & tissue.index.isin(expression.index)]
        if len(samples) > 0:
            blocks[cohort.id] = expression.loc[samples]
    genes = common_genes(blocks)
    centered = {cid: b[genes] - b[genes].mean() for cid, b in blocks.items()}
    return pd.concat(centered, names=["cohort", "sample"])


def _select_and_scale(train: pd.DataFrame, other: pd.DataFrame | None, n_genes: int):
    """Gènes les plus variables et mise à l'échelle, calculés sur `train` seulement."""
    genes = train.var().sort_values(ascending=False).index[:n_genes]
    mean, std = train[genes].mean(), train[genes].std(ddof=0).replace(0, 1.0)
    scaled_train = (train[genes] - mean) / std
    scaled_other = None if other is None else (other[genes] - mean) / std
    return scaled_train, scaled_other


def run_subtypes(
    data: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    cohorts: list[Cohort],
    responses: dict[str, pd.Series],
    k_values: tuple[int, ...] = (2, 3, 4, 5),
    n_genes: int = 2000,
    n_resamples: int = 100,
    log: Callable[[str], None] | None = None,
) -> dict[str, pd.DataFrame]:
    """Sous-types des lésions, leur stabilité entre cohortes et leur lien avec la réponse."""
    from scipy.stats import fisher_exact
    from sklearn.metrics import adjusted_rand_score

    from hs_stratify.subtypes import (
        align_labels,
        best_k,
        consensus_clustering,
        marker_genes,
        nearest_centroid,
    )

    say = log or (lambda m: None)
    lesions = lesional_matrix(data, cohorts)
    say(f"{lesions.shape[0]} lésions, {lesions.shape[1]} gènes communs")

    X, _ = _select_and_scale(lesions, None, n_genes)
    say("Clustering de consensus sur toutes les lésions...")
    results = consensus_clustering(X.to_numpy(), k_values, n_resamples=n_resamples)
    k = best_k(results)
    labels = pd.Series(results[k].labels, index=X.index, name="subtype")
    pac = pd.DataFrame({"k": list(results), "pac": [round(r.pac, 3) for r in results.values()]})
    say(f"k retenu : {k}")

    # Stabilité : les sous-types trouvés sans une cohorte se retrouvent-ils sur cette cohorte ?
    stability, predicted_all = [], []
    groups = lesions.index.get_level_values("cohort")
    for held_out in dict.fromkeys(groups):
        say(f"  stabilité : sous-types appris sans {held_out}")
        train_raw, test_raw = lesions[groups != held_out], lesions[groups == held_out]
        train, test = _select_and_scale(train_raw, test_raw, n_genes)
        fold = consensus_clustering(train.to_numpy(), (k,), n_resamples=n_resamples // 2)[k]
        # Numéros alignés sur les sous-types globaux via les seuls échantillons d'entraînement.
        fold_labels = align_labels(fold.labels, labels.loc[train.index].to_numpy())
        predicted = nearest_centroid(train.to_numpy(), fold_labels, test.to_numpy())
        reference = labels.loc[test.index].to_numpy()
        predicted_all.append(pd.Series(predicted, index=test.index))
        stability.append(
            {
                "cohort": held_out,
                "n_lesions": len(test),
                "ari": round(adjusted_rand_score(reference, predicted), 3),
                "agreement": round(float(np.mean(reference == predicted)), 3),
            }
        )
    pooled = pd.concat(predicted_all).loc[labels.index]
    stability.append(
        {
            "cohort": "all",
            "n_lesions": len(labels),
            "ari": round(adjusted_rand_score(labels, pooled), 3),
            "agreement": round(float(np.mean(labels.to_numpy() == pooled.to_numpy())), 3),
        }
    )

    composition = pd.crosstab(groups, labels.to_numpy(), rownames=["cohort"], colnames=["subtype"])

    # Lien avec la réponse au traitement, sur les seuls patients étiquetés.
    rows = []
    for cohort_id, response in responses.items():
        for sample, value in response.items():
            if (cohort_id, sample) in labels.index:
                rows.append(
                    {
                        "cohort": cohort_id,
                        "sample": sample,
                        "subtype": labels.loc[(cohort_id, sample)],
                        "response": int(value),
                    }
                )
    response_table = pd.DataFrame(rows, columns=["cohort", "sample", "subtype", "response"])
    response_summary = pd.DataFrame()
    if not response_table.empty:
        response_summary = (
            response_table.groupby("subtype")["response"]
            .agg(patients="size", responders="sum")
            .assign(rate=lambda t: (t["responders"] / t["patients"]).round(3))
        )
        if k == 2 and len(response_summary) == 2:
            table = pd.crosstab(response_table["subtype"], response_table["response"])
            table = table.reindex(index=[0, 1], columns=[0, 1], fill_value=0)
            _, p_value = fisher_exact(table.to_numpy())
            response_summary["fisher_p"] = round(float(p_value), 4)

    return {
        "pac": pac,
        "labels": labels.reset_index(),
        "composition": composition,
        "stability": pd.DataFrame(stability),
        "response": response_table,
        "response_summary": response_summary,
        "markers": marker_genes(X, labels.to_numpy()),
    }


def response_dataset(
    data: dict[str, tuple[pd.DataFrame, pd.DataFrame]], responses: dict[str, pd.Series]
) -> tuple[pd.DataFrame, pd.Series]:
    """Biopsies étiquetées, gènes communs, chaque échantillon centré-réduit sur ses gènes."""
    matrices, targets = {}, {}
    for cohort_id, response in responses.items():
        if cohort_id not in data:
            continue
        expression, _ = data[cohort_id]
        response = response[response.index.isin(expression.index)]
        if response.empty:
            continue
        matrices[cohort_id] = expression.loc[response.index]
        targets[cohort_id] = response.astype(int)
    genes = common_genes(matrices)
    X = pd.concat(
        {cid: standardize_per_sample(m[genes]) for cid, m in matrices.items()},
        names=["cohort", "sample"],
    )
    y = pd.concat(targets, names=["cohort", "sample"]).loc[X.index]
    return X, y


def run_response(
    data: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    responses: dict[str, pd.Series],
    n_genes: int = 2000,
    C: float = 0.1,
    n_permutations: int = 200,
    n_bootstrap: int = 1000,
    n_loo_permutations: int = 50,
    log: Callable[[str], None] | None = None,
) -> pd.DataFrame:
    """Prédiction de la réponse : d'une cohorte vers l'autre, puis en leave-one-out groupé.

    Le leave-one-out mélange les deux cohortes : il est donné pour comparaison, il est
    optimiste puisque le modèle voit des patients de la même étude que celui qu'il prédit.
    """
    from sklearn.metrics import roc_auc_score

    from hs_stratify.response import (
        bootstrap_auroc_ci,
        fit_score,
        leave_one_out_p_value,
        leave_one_out_scores,
        permutation_p_value,
    )

    say = log or (lambda m: None)
    X, y = response_dataset(data, responses)
    say(f"{len(y)} patients étiquetés, {X.shape[1]} gènes communs")

    def make_model():
        return ridge_model(C=C, n_genes=n_genes)

    def counts(mask) -> str:
        values = y[mask]
        return f"{int(values.sum())}R/{int((values == 0).sum())}NR"

    groups = X.index.get_level_values("cohort")
    rows = []
    for train_id in dict.fromkeys(groups):
        for test_id in dict.fromkeys(groups):
            if train_id == test_id:
                continue
            say(f"  entraînement sur {train_id}, test sur {test_id}...")
            train, test = groups == train_id, groups == test_id
            X_tr, y_tr = X[train].to_numpy(), y[train].to_numpy()
            X_te, y_te = X[test].to_numpy(), y[test].to_numpy()
            scores = fit_score(make_model, X_tr, y_tr, X_te)
            auroc = roc_auc_score(y_te, scores)
            low, high = bootstrap_auroc_ci(y_te, scores, n_bootstrap)
            p = permutation_p_value(make_model, X_tr, y_tr, X_te, y_te, auroc, n_permutations)
            rows.append(
                {
                    "setting": "cross-cohort",
                    "train": train_id,
                    "test": test_id,
                    "n_train": counts(train),
                    "n_test": counts(test),
                    "auroc": auroc,
                    "ci_low": low,
                    "ci_high": high,
                    "perm_p": p,
                }
            )

    say("  leave-one-out sur les deux cohortes réunies (optimiste)...")
    X_all, y_all = X.to_numpy(), y.to_numpy()
    scores = leave_one_out_scores(make_model, X_all, y_all)
    auroc = roc_auc_score(y_all, scores)
    low, high = bootstrap_auroc_ci(y_all, scores, n_bootstrap)
    p = leave_one_out_p_value(make_model, X_all, y_all, auroc, n_loo_permutations)
    everyone = np.ones(len(y), dtype=bool)
    rows.append(
        {
            "setting": "pooled leave-one-out",
            "train": "all but one",
            "test": "each patient",
            "n_train": counts(everyone),
            "n_test": counts(everyone),
            "auroc": auroc,
            "ci_low": low,
            "ci_high": high,
            "perm_p": p,
        }
    )
    return pd.DataFrame(rows).round(3)
