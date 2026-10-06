"""Registre des cohortes GEO utilisées dans le projet."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "configs" / "cohorts.yaml"
ROLES = {"discovery", "validation"}
PLATFORMS = {"rnaseq", "array"}


@dataclass(frozen=True)
class Cohort:
    id: str
    role: str
    platform: str
    response_labels: str = "none"
    # Règles pour extraire la réponse au traitement des métadonnées (voir labels.py).
    labels: dict | None = None
    # Règles pour lire le type de tissu (lésion, peau saine...) dans les métadonnées.
    tissue: dict | None = None


def load_cohorts(path: Path | str = DEFAULT_CONFIG) -> list[Cohort]:
    """Lit le fichier de configuration et renvoie la liste des cohortes."""
    with open(path, encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    cohorts = [Cohort(**entry) for entry in raw["cohorts"]]
    for cohort in cohorts:
        if not cohort.id.startswith("GSE"):
            raise ValueError(f"Identifiant GEO invalide : {cohort.id}")
        if cohort.role not in ROLES:
            raise ValueError(f"Rôle inconnu pour {cohort.id} : {cohort.role}")
        if cohort.platform not in PLATFORMS:
            raise ValueError(f"Plateforme inconnue pour {cohort.id} : {cohort.platform}")
        if cohort.labels is not None and cohort.response_labels == "none":
            raise ValueError(f"{cohort.id} : règles d'étiquetage sans response_labels")

    ids = [c.id for c in cohorts]
    if len(ids) != len(set(ids)):
        raise ValueError("Une cohorte apparaît plusieurs fois dans la configuration")
    return cohorts
