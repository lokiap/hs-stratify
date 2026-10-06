import pytest

from hs_stratify.cohorts import load_cohorts


def test_default_config_loads():
    cohorts = load_cohorts()
    roles = {c.role for c in cohorts}
    assert roles == {"discovery", "validation"}
    assert len(cohorts) == 9
    assert {c.platform for c in cohorts} == {"rnaseq", "array"}
    labelled = [c.id for c in cohorts if c.response_labels != "none"]
    assert labelled == ["GSE155176", "GSE213761"]


def test_rejects_unknown_role(tmp_path):
    config = tmp_path / "cohorts.yaml"
    config.write_text(
        "cohorts:\n  - id: GSE1\n    role: other\n    platform: rnaseq\n", encoding="utf-8"
    )
    with pytest.raises(ValueError):
        load_cohorts(config)


def test_rejects_duplicates(tmp_path):
    config = tmp_path / "cohorts.yaml"
    config.write_text(
        "cohorts:\n  - id: GSE1\n    role: discovery\n    platform: rnaseq\n"
        "  - id: GSE1\n    role: validation\n    platform: array\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        load_cohorts(config)
