import pytest

from hs_stratify.cohorts import load_cohorts


def test_default_config_loads():
    cohorts = load_cohorts()
    roles = {c.role for c in cohorts}
    assert roles == {"discovery", "validation"}
    assert len(cohorts) == 9


def test_rejects_unknown_role(tmp_path):
    config = tmp_path / "cohorts.yaml"
    config.write_text("cohorts:\n  - id: GSE1\n    role: other\n", encoding="utf-8")
    with pytest.raises(ValueError):
        load_cohorts(config)


def test_rejects_duplicates(tmp_path):
    config = tmp_path / "cohorts.yaml"
    config.write_text(
        "cohorts:\n  - id: GSE1\n    role: discovery\n  - id: GSE1\n    role: validation\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        load_cohorts(config)
