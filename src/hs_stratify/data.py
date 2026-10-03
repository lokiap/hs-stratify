"""Téléchargement des cohortes depuis GEO.

Nécessite l'extra `data` : pip install -e ".[data]"
"""

from __future__ import annotations

from pathlib import Path

RAW_DIR = Path("data/raw")


def download_cohort(gse_id: str, dest: Path = RAW_DIR):
    """Télécharge une série GEO (SOFT) dans `dest` et renvoie l'objet GEOparse."""
    try:
        import GEOparse
    except ImportError as exc:
        raise ImportError("GEOparse n'est pas installé : pip install -e '.[data]'") from exc

    dest.mkdir(parents=True, exist_ok=True)
    return GEOparse.get_GEO(geo=gse_id, destdir=str(dest), silent=True)


def build_expression_matrix(gse_id: str, dest: Path = RAW_DIR):
    """Construit la matrice gènes x échantillons d'une cohorte.

    À faire : gérer séparément puces (GPL Affymetrix/Illumina) et RNA-seq,
    puis ramener les sondes au niveau du gène.
    """
    raise NotImplementedError
