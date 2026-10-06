import pandas as pd
import pytest

from hs_stratify.cohorts import load_cohorts
from hs_stratify.labels import cohort_labels as response_labels_for
from hs_stratify.labels import response_labels

META = pd.DataFrame(
    {
        "ch1:tissue type": ["Lesional", "Lesional", "Lesional", "Lesional", "Non-lesional"],
        "ch1:treatment category": ["PRE-TNF-alpha inhibitor"] * 3 + ["ON-TNF-alpha inhibitor"] * 2,
        "ch1:response": ["Yes", "no", "Partial", "Yes", "yes"],
    },
    index=["GSM1", "GSM2", "GSM3", "GSM4", "GSM5"],
)


def test_gse155176_rules_keep_baseline_lesions_with_clear_response():
    rules = next(c.labels for c in load_cohorts() if c.id == "GSE155176")
    labels = response_labels(META, rules)
    assert labels.to_dict() == {"GSM1": 1, "GSM2": 0}


def test_missing_field_is_reported():
    rules = {"response": {"field": "ch1:hiscr", "positive": ["yes"], "negative": ["no"]}}
    with pytest.raises(KeyError):
        response_labels(META, rules)


def test_gse213761_hiscr_from_lesion_counts():
    rules = next(c.labels for c in load_cohorts() if c.id == "GSE213761")
    # Patient 1 : AN 6 -> 2, répondeur. Patient 2 : AN 4 -> 3, non-répondeur.
    # Patient 3 : AN 6 -> 1 mais fistules en hausse, non-répondeur.
    # Patient 4 : AN 2 au départ, écarté. Patient 5 : pas de visite V12, écarté.
    rows = []
    for subject, v1, v12 in [
        ("1", (2, 4, 1), (0, 2, 1)),
        ("2", (1, 3, 0), (1, 2, 0)),
        ("3", (3, 3, 1), (0, 1, 2)),
        ("4", (0, 2, 0), (0, 0, 0)),
        ("5", (3, 3, 0), None),
    ]:
        for visit, counts in [("V1", v1), ("V12", v12)]:
            if counts is None:
                continue
            for site in ["LS", "NL"]:
                abscess, nodule, fistula = counts
                rows.append(
                    {
                        "ch1:subject": subject,
                        "ch1:visit": visit,
                        "ch1:site": site,
                        "ch1:abcess count": str(abscess),
                        "ch1:node count": str(nodule),
                        "ch1:fistula count": str(fistula),
                    }
                )
    meta = pd.DataFrame(rows, index=[f"GSM{i}" for i in range(len(rows))])

    labels = response_labels_for(meta, rules)

    by_subject = meta.loc[labels.index, "ch1:subject"].map(str).tolist()
    assert dict(zip(by_subject, labels, strict=True)) == {"1": 1, "2": 0, "3": 0}
    assert set(meta.loc[labels.index, "ch1:site"]) == {"LS"}
    assert set(meta.loc[labels.index, "ch1:visit"]) == {"V1"}
