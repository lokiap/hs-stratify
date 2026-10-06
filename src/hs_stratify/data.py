"""Construction des matrices d'expression de chaque cohorte à partir de GEO."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from hs_stratify import geo
from hs_stratify.cohorts import Cohort
from hs_stratify.preprocess import ensure_log2, log_cpm, probes_to_genes

RAW_DIR = Path("data/raw")


def fetch_metadata(gse_id: str, raw_dir: Path = RAW_DIR):
    """Télécharge le series matrix et renvoie (métadonnées, expression ou None, plateformes)."""
    path = geo.download(geo.series_matrix_url(gse_id), raw_dir / f"{gse_id}_series_matrix.txt.gz")
    return geo.parse_series_matrix(path)


def load_cohort(cohort: Cohort, raw_dir: Path = RAW_DIR) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Renvoie (expression log2 échantillons x gènes, métadonnées) pour une cohorte."""
    meta, array_values, platforms = fetch_metadata(cohort.id, raw_dir)

    if cohort.platform == "rnaseq":
        counts_path = geo.download(
            geo.ncbi_counts_url(cohort.id), raw_dir / f"{cohort.id}_ncbi_counts.tsv.gz"
        )
        annot_path = geo.download(geo.ncbi_annot_url(), raw_dir / geo.HUMAN_ANNOT)
        genes_x_samples = log_cpm(geo.read_ncbi_counts(counts_path, annot_path))
    else:
        if array_values is None:
            raise ValueError(f"{cohort.id} : pas de valeurs d'expression dans le series matrix")
        if len(platforms) != 1:
            raise ValueError(f"{cohort.id} : plusieurs plateformes {platforms}, cas non géré")
        annot_path = geo.download(
            geo.platform_annot_url(platforms[0]), raw_dir / f"{platforms[0]}.annot.gz"
        )
        genes_x_samples = probes_to_genes(
            ensure_log2(array_values), geo.read_platform_annot(annot_path)
        )

    samples = genes_x_samples.columns.intersection(meta.index)
    return genes_x_samples[samples].T, meta.loc[samples]
