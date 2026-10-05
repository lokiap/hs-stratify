"""Découpages d'évaluation : on teste toujours sur une cohorte jamais vue."""

from __future__ import annotations

from collections.abc import Iterator, Sequence

import numpy as np


def leave_one_cohort_out(
    cohorts: Sequence[str],
) -> Iterator[tuple[str, np.ndarray, np.ndarray]]:
    """Pour chaque cohorte, renvoie (cohorte testée, indices d'entraînement, indices de test).

    `cohorts` donne la cohorte de chaque échantillon, dans l'ordre des lignes.
    """
    labels = np.asarray(cohorts)
    unique = list(dict.fromkeys(labels.tolist()))
    if len(unique) < 2:
        raise ValueError("Il faut au moins deux cohortes pour un découpage par cohorte")

    for held_out in unique:
        test = np.flatnonzero(labels == held_out)
        train = np.flatnonzero(labels != held_out)
        yield held_out, train, test
