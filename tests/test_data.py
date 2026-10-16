"""Chargement d'une cohorte de bout en bout, hors ligne : les fichiers sont déjà en cache."""

import gzip

from hs_stratify import geo
from hs_stratify.cohorts import Cohort
from hs_stratify.data import load_cohort

from .test_geo import SERIES_MATRIX


def _write_gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)


def test_load_array_cohort(tmp_path):
    _write_gz(tmp_path / "GSE1_series_matrix.txt.gz", SERIES_MATRIX)
    _write_gz(tmp_path / "GPL570.annot.gz", "ID\tGene symbol\np1\tTNF\np2\tKRT16\n")

    expression, meta = load_cohort(Cohort("GSE1", "validation", "array"), raw_dir=tmp_path)

    assert expression.shape == (3, 2)
    assert set(expression.columns) == {"TNF", "KRT16"}
    assert list(meta.index) == list(expression.index)


def test_load_rnaseq_cohort(tmp_path):
    matrix = SERIES_MATRIX.split("!series_matrix_table_begin")[0]
    _write_gz(tmp_path / "GSE1_series_matrix.txt.gz", matrix)
    _write_gz(
        tmp_path / "GSE1_ncbi_counts.tsv.gz",
        "GeneID\tGSM1\tGSM2\tGSM3\n7124\t500\t800\t300\n3605\t900\t100\t700\n",
    )
    _write_gz(tmp_path / geo.HUMAN_ANNOT, "GeneID\tSymbol\n7124\tTNF\n3605\tIL17A\n")

    expression, meta = load_cohort(Cohort("GSE1", "discovery", "rnaseq"), raw_dir=tmp_path)

    assert expression.shape == (3, 2)
    assert meta.loc["GSM1", "ch1:response"] == "responder"
