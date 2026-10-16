import gzip
import io

import pandas as pd
import pytest

from hs_stratify import geo

SERIES_MATRIX = """!Series_title\t"HS skin"
!Series_platform_id\t"GPL570"
!Sample_title\t"HS_1"\t"HS_2"\t"Ctrl_1"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_source_name_ch1\t"skin"\t"skin"\t"skin"
!Sample_characteristics_ch1\t"tissue: lesional"\t"tissue: lesional"\t"tissue: healthy"
!Sample_characteristics_ch1\t"response: responder"\t"response: non-responder"\t""
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"
"p1"\t7.5\t8.1\t6.0
"p2"\t5.0\t5.5\t4.9
!series_matrix_table_end
"""


def test_stub():
    assert geo._stub("GSE154773") == "GSE154nnn"
    assert geo._stub("GSE1234") == "GSE1nnn"
    assert geo._stub("GSE123") == "GSEnnn"
    assert geo.series_matrix_url("GSE72702").endswith(
        "/series/GSE72nnn/GSE72702/matrix/GSE72702_series_matrix.txt.gz"
    )


def test_parse_series_matrix(tmp_path):
    path = tmp_path / "GSE1_series_matrix.txt.gz"
    with gzip.open(path, "wt") as f:
        f.write(SERIES_MATRIX)

    meta, expression, platforms = geo.parse_series_matrix(path)

    assert platforms == ["GPL570"]
    assert list(meta.index) == ["GSM1", "GSM2", "GSM3"]
    assert list(meta["ch1:tissue"]) == ["lesional", "lesional", "healthy"]
    assert meta.loc["GSM2", "ch1:response"] == "non-responder"
    assert pd.isna(meta.loc["GSM3", "ch1:response"])
    assert expression.shape == (2, 3)
    assert expression.loc["p1", "GSM2"] == 8.1


def test_read_platform_annot_drops_ambiguous_probes(tmp_path):
    path = tmp_path / "GPL1.annot"
    path.write_text(
        "!Annotation_platform = GPL1\n"
        "ID\tGene symbol\tGene ID\n"
        "p1\tTNF\t7124\n"
        "p2\tIL17A///IL17F\t3605\n"
        "p3\t\t\n"
    )
    symbols = geo.read_platform_annot(path)
    assert symbols.to_dict() == {"p1": "TNF"}


def test_read_ncbi_counts_sums_by_symbol(tmp_path):
    counts = tmp_path / "counts.tsv"
    counts.write_text("GeneID\tGSM1\tGSM2\n7124\t10\t20\n3605\t1\t2\n9999\t5\t5\n")
    # En-tête du fichier gene_info du FTP NCBI : la première colonne commence par « # ».
    annot = tmp_path / "annot.tsv"
    annot.write_text(
        "#tax_id\tGeneID\tSymbol\tLocusTag\n"
        "9606\t7124\tTNF\t-\n"
        "9606\t3605\tIL17A\t-\n"
        "9606\t1\tA1BG\t-\n"
    )

    result = geo.read_ncbi_counts(counts, annot)

    assert sorted(result.index) == ["IL17A", "TNF"]
    assert result.loc["TNF", "GSM2"] == 20


def test_download_replaces_cached_html_and_reports_error(tmp_path, monkeypatch):
    dest = tmp_path / "annot.tsv.gz"
    dest.write_text("<!DOCTYPE html><html><title>Error</title></html>")

    def fake_urlopen(request, timeout):
        return io.BytesIO(b"<!DOCTYPE html><html><title>Too Many Requests</title></html>")

    monkeypatch.setattr(geo.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(geo.time, "sleep", lambda s: None)

    with pytest.raises(RuntimeError, match="Too Many Requests"):
        geo.download("https://example.org/annot", dest)
    assert not dest.exists()


def test_download_keeps_valid_cache(tmp_path, monkeypatch):
    dest = tmp_path / "ok.tsv.gz"
    with gzip.open(dest, "wt") as f:
        f.write("GeneID\tSymbol\n")
    monkeypatch.setattr(geo.urllib.request, "urlopen", lambda *a, **k: 1 / 0)
    assert geo.download("https://example.org/ok", dest) == dest


def test_check_gzip_names_recaptcha(tmp_path):
    path = tmp_path / "x.gz"
    path.write_text('<!doctype html><html><head><base href="https://www.google.com/recaptcha/">')
    with pytest.raises(ValueError, match="anti-robot"):
        geo._check_gzip(path)


def test_counts_urls_prefer_ftp():
    first, second = geo.ncbi_counts_urls("GSE151243")
    assert first == (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151243/suppl/"
        "GSE151243_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    )
    assert second.startswith("https://www.ncbi.nlm.nih.gov/geo/download/")


def test_download_first_falls_back(tmp_path, monkeypatch):
    dest = tmp_path / "counts.tsv.gz"
    calls = []

    def fake_urlopen(request, timeout):
        calls.append(request.full_url)
        if request.full_url.startswith("https://a/"):
            raise OSError("404")
        buffer = io.BytesIO()
        with gzip.GzipFile(fileobj=buffer, mode="wb") as f:
            f.write(b"GeneID\tGSM1\n")
        return io.BytesIO(buffer.getvalue())

    monkeypatch.setattr(geo.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(geo.time, "sleep", lambda s: None)

    assert geo.download_first(["https://a/x", "https://b/x"], dest) == dest
    assert calls[-1] == "https://b/x"
