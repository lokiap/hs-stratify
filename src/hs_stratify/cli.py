"""Interface en ligne de commande."""

from __future__ import annotations

import argparse

from hs_stratify.cohorts import load_cohorts
from hs_stratify.data import download_cohort


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="hs-stratify")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("cohorts", help="Lister les cohortes configurées")
    dl = sub.add_parser("download", help="Télécharger les cohortes depuis GEO")
    dl.add_argument("ids", nargs="*", help="Identifiants GSE (toutes par défaut)")
    args = parser.parse_args(argv)

    cohorts = load_cohorts()
    if args.command == "cohorts":
        for c in cohorts:
            print(f"{c.id}\t{c.role}\tréponse: {c.response_labels}")
    elif args.command == "download":
        for gse_id in args.ids or [c.id for c in cohorts]:
            print(f"Téléchargement de {gse_id}...")
            download_cohort(gse_id)


if __name__ == "__main__":
    main()
