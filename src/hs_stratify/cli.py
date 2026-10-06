"""Interface en ligne de commande."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from hs_stratify.cohorts import load_cohorts
from hs_stratify.data import fetch_metadata, load_cohort
from hs_stratify.experiments import load_processed, run_lesional, summarize
from hs_stratify.labels import cohort_labels

PROCESSED_DIR = Path("data/processed")
RESULTS_DIR = Path("results")


def _inspect(gse_id: str) -> None:
    """Affiche les champs de métadonnées d'une série et leurs valeurs les plus fréquentes."""
    meta, _, platforms = fetch_metadata(gse_id)
    print(f"{gse_id} : {len(meta)} échantillons, plateforme(s) {', '.join(platforms)}")
    for column in meta.columns:
        if column.startswith("ch1:") or column in {"title", "source_name_ch1"}:
            counts = meta[column].value_counts()
            preview = ", ".join(f"{v} ({n})" for v, n in counts.head(6).items())
            print(f"  {column} [{counts.size} valeurs] : {preview}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="hs-stratify")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("cohorts", help="Lister les cohortes configurées")
    inspect = sub.add_parser("inspect", help="Afficher les métadonnées d'une série GEO")
    inspect.add_argument("id", help="Identifiant GSE")
    build = sub.add_parser("build", help="Télécharger et construire les matrices d'expression")
    build.add_argument("ids", nargs="*", help="Identifiants GSE (toutes par défaut)")
    evaluate = sub.add_parser("evaluate", help="Lancer une expérience sur les données construites")
    evaluate.add_argument("experiment", choices=["lesional"], help="Expérience à lancer")
    args = parser.parse_args(argv)

    cohorts = load_cohorts()
    if args.command == "cohorts":
        for c in cohorts:
            print(f"{c.id}\t{c.role}\t{c.platform}\tréponse: {c.response_labels}")
    elif args.command == "inspect":
        _inspect(args.id)
    elif args.command == "build":
        wanted = set(args.ids) or {c.id for c in cohorts}
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        for cohort in (c for c in cohorts if c.id in wanted):
            print(f"{cohort.id} : construction...")
            expression, meta = load_cohort(cohort)
            expression.to_csv(PROCESSED_DIR / f"{cohort.id}_expression.csv.gz")
            meta.to_csv(PROCESSED_DIR / f"{cohort.id}_metadata.csv.gz")
            print(f"  {expression.shape[0]} échantillons x {expression.shape[1]} gènes")
            if cohort.labels:
                labels = cohort_labels(meta, cohort.labels)
                labels.to_csv(PROCESSED_DIR / f"{cohort.id}_labels.csv")
                n_pos = int(labels.sum())
                print(f"  réponse : {n_pos} répondeurs, {len(labels) - n_pos} non-répondeurs")
    elif args.command == "evaluate":
        data = {
            c.id: load_processed(c.id, PROCESSED_DIR)
            for c in cohorts
            if (PROCESSED_DIR / f"{c.id}_expression.csv.gz").exists()
        }
        if not data:
            parser.error("aucune donnée construite : lancer d'abord `hs-stratify build`")
        results = run_lesional(data, cohorts)
        RESULTS_DIR.mkdir(exist_ok=True)
        results.to_csv(RESULTS_DIR / "lesional_loco.csv", index=False)
        with pd.option_context("display.width", 120):
            print(results.round(3).to_string(index=False))
            print()
            print(summarize(results).to_string())


if __name__ == "__main__":
    main()
