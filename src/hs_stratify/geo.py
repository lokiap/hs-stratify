"""Accès aux fichiers GEO (NCBI) et lecture de leurs formats.

Trois sources sont utilisées, toutes publiques et sans compte :
- le « series matrix » de chaque série : métadonnées des échantillons, et pour les puces,
  les valeurs d'expression par sonde ;
- les comptages RNA-seq recalculés par NCBI sur GRCh38, identiques d'une étude à l'autre
  (plus simple et plus homogène que les fichiers déposés par chaque équipe) ;
- l'annotation des plateformes de puces (sonde -> symbole de gène).
"""

from __future__ import annotations

import gzip
import io
import shutil
import urllib.request
from pathlib import Path

import pandas as pd

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
GEO_DOWNLOAD = "https://www.ncbi.nlm.nih.gov/geo/download/"
HUMAN_ANNOT = "Human.GRCh38.p13.annot.tsv.gz"


def _stub(accession: str) -> str:
    """GSE154773 -> GSE154nnn, GSE1234 -> GSE1nnn, GSE123 -> GSEnnn."""
    prefix = accession[:3]
    digits = accession[3:]
    return prefix + (digits[:-3] if len(digits) > 3 else "") + "nnn"


def series_matrix_url(gse_id: str) -> str:
    return f"{GEO_FTP}/series/{_stub(gse_id)}/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"


def platform_annot_url(gpl_id: str) -> str:
    return f"{GEO_FTP}/platforms/{_stub(gpl_id)}/{gpl_id}/annot/{gpl_id}.annot.gz"


def ncbi_counts_url(gse_id: str) -> str:
    file = f"{gse_id}_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    return f"{GEO_DOWNLOAD}?type=rnaseq_counts&acc={gse_id}&format=file&file={file}"


def ncbi_annot_url() -> str:
    return f"{GEO_DOWNLOAD}?format=file&type=rnaseq_counts&file={HUMAN_ANNOT}"


def download(url: str, dest: Path) -> Path:
    """Télécharge `url` vers `dest`, sauf si le fichier est déjà en cache."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "hs-stratify"})
    with urllib.request.urlopen(request, timeout=120) as response, open(tmp, "wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(dest)
    return dest


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _split_fields(line: str) -> list[str]:
    return [field.strip().strip('"') for field in line.rstrip("\n").split("\t")]


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame | None, list[str]]:
    """Lit un series matrix GEO.

    Renvoie (métadonnées échantillons x champs, expression sondes x échantillons ou None,
    identifiants de plateforme). Les lignes « !Sample_characteristics_ch1 » de la forme
    « clé: valeur » deviennent une colonne par clé.
    """
    sample_fields: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    platforms: list[str] = []
    table_lines: list[str] = []
    in_table = False

    with _open_text(path) as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
                continue
            if line.startswith("!Series_platform_id"):
                platforms.extend(_split_fields(line)[1:])
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(_split_fields(line)[1:])
            elif line.startswith("!Sample_"):
                fields = _split_fields(line)
                key = fields[0][len("!Sample_") :]
                if key not in sample_fields:
                    sample_fields[key] = fields[1:]

    accessions = sample_fields.get("geo_accession")
    if not accessions:
        raise ValueError(f"Aucun échantillon trouvé dans {path}")

    meta = pd.DataFrame(
        {k: v for k, v in sample_fields.items() if len(v) == len(accessions)},
        index=pd.Index(accessions, name="sample"),
    )
    for row in characteristics:
        for sample, cell in zip(accessions, row, strict=False):
            if ":" in cell:
                key, value = cell.split(":", 1)
                meta.loc[sample, f"ch1:{key.strip().lower()}"] = value.strip()

    expression = None
    if len(table_lines) > 1:
        expression = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
        expression.index = expression.index.astype(str)
        expression = expression.apply(pd.to_numeric, errors="coerce")
    return meta, expression, platforms


def read_platform_annot(path: Path) -> pd.Series:
    """Lit un fichier .annot GEO et renvoie la correspondance sonde -> symbole de gène."""
    with _open_text(path) as f:
        lines = [line for line in f if not line.startswith(("!", "#", "^"))]
    table = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    symbols = table.set_index("ID")["Gene symbol"].dropna()
    # Une sonde ambiguë (« A///B ») est écartée plutôt qu'attribuée au hasard.
    return symbols[~symbols.str.contains("///", regex=False)]


def read_ncbi_counts(counts_path: Path, annot_path: Path) -> pd.DataFrame:
    """Lit les comptages NCBI (GeneID x échantillons) et indexe par symbole de gène."""
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    annot = pd.read_csv(annot_path, sep="\t", index_col=0, usecols=["GeneID", "Symbol"])
    counts = counts.join(annot, how="inner").dropna(subset=["Symbol"])
    return counts.groupby("Symbol").sum()
