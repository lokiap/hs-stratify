"""Extraction de la réponse au traitement à partir des métadonnées GEO."""

from __future__ import annotations

import pandas as pd


def _normalize(values: pd.Series) -> pd.Series:
    return values.astype("string").str.strip().str.lower()


def response_labels(meta: pd.DataFrame, rules: dict) -> pd.Series:
    """Applique les règles d'une cohorte et renvoie 1 (répondeur) / 0 (non-répondeur).

    Seuls les échantillons qui passent tous les filtres et dont la réponse est explicitement
    positive ou négative sont gardés : les valeurs ambiguës (NA, Partial...) sont écartées.
    """
    keep = pd.Series(True, index=meta.index)
    for field, allowed in rules.get("filters", {}).items():
        if field not in meta.columns:
            raise KeyError(f"Champ de filtre absent des métadonnées : {field}")
        keep &= _normalize(meta[field]).isin([str(v).lower() for v in allowed]).fillna(False)

    response = rules["response"]
    if response["field"] not in meta.columns:
        raise KeyError(f"Champ de réponse absent des métadonnées : {response['field']}")
    values = _normalize(meta.loc[keep, response["field"]])
    positive = values.isin([str(v).lower() for v in response["positive"]]).fillna(False)
    negative = values.isin([str(v).lower() for v in response["negative"]]).fillna(False)

    labels = pd.Series(pd.NA, index=values.index, dtype="Int64")
    labels[positive] = 1
    labels[negative] = 0
    return labels.dropna().astype(int).rename("response")


def hiscr_labels(meta: pd.DataFrame, rules: dict) -> pd.Series:
    """Calcule la réponse HiSCR à partir des comptes de lésions avant et après traitement.

    HiSCR (Kimball et al., 2014) : répondeur si le nombre d'abcès + nodules inflammatoires (AN)
    baisse d'au moins 50 %, sans hausse du nombre d'abcès ni de fistules qui suintent.
    Les comptes sont propres au patient et à la visite : on les lit sur tous ses échantillons,
    puis on attribue l'étiquette aux échantillons retenus par les filtres (biopsies initiales).
    Les patients avec moins de `min_baseline_an` lésions au départ sont écartés, comme dans
    les essais cliniques.
    """
    h = rules["hiscr"]
    counts = meta[[h["subject"], h["visit"], h["abscess"], h["nodule"], h["fistula"]]].copy()
    counts.columns = ["subject", "visit", "abscess", "nodule", "fistula"]
    for col in ["abscess", "nodule", "fistula"]:
        counts[col] = pd.to_numeric(counts[col], errors="coerce")
    per_visit = counts.dropna().groupby(["subject", "visit"]).first()

    visits = set(per_visit.index.get_level_values("visit"))
    for name in (h["baseline"], h["followup"]):
        if name not in visits:
            raise ValueError(f"Visite {name} absente des métadonnées")
    base = per_visit.xs(h["baseline"], level="visit")
    follow = per_visit.xs(h["followup"], level="visit")
    both = base.join(follow, lsuffix="_0", rsuffix="_1", how="inner")
    an_0 = both["abscess_0"] + both["nodule_0"]
    an_1 = both["abscess_1"] + both["nodule_1"]
    both = both[an_0 >= h.get("min_baseline_an", 3)]
    an_0, an_1 = an_0[both.index], an_1[both.index]

    responder = (
        (an_1 <= 0.5 * an_0)
        & (both["abscess_1"] <= both["abscess_0"])
        & (both["fistula_1"] <= both["fistula_0"])
    ).astype(int)

    keep = pd.Series(True, index=meta.index)
    for field, allowed in rules.get("filters", {}).items():
        keep &= _normalize(meta[field]).isin([str(v).lower() for v in allowed]).fillna(False)
    subjects = meta.loc[keep, h["subject"]]
    labels = subjects.map(responder).dropna().astype(int)
    return labels.rename("response")


def cohort_labels(meta: pd.DataFrame, rules: dict) -> pd.Series:
    """Choisit la bonne méthode selon les règles de la cohorte."""
    return hiscr_labels(meta, rules) if "hiscr" in rules else response_labels(meta, rules)
