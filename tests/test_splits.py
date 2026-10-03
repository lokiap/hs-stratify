import numpy as np
import pytest

from hs_stratify.splits import leave_one_cohort_out


def test_each_cohort_held_out_once_without_leakage():
    cohorts = ["A", "A", "B", "C", "C", "C"]
    folds = list(leave_one_cohort_out(cohorts))

    assert [held for held, _, _ in folds] == ["A", "B", "C"]
    for held, train, test in folds:
        assert set(np.asarray(cohorts)[test]) == {held}
        assert held not in set(np.asarray(cohorts)[train])
        assert len(train) + len(test) == len(cohorts)


def test_needs_two_cohorts():
    with pytest.raises(ValueError):
        list(leave_one_cohort_out(["A", "A"]))
