# hs-stratify

Predicting adalimumab response in HS (Hidradenitis Suppurativa) from skin transcriptomics, with cross-cohort validation.

[English](#english) · [Français](#français)

> **Status:** work in progress. The repository structure, cohort registry and evaluation protocol are in place; data loading and models are not yet implemented.

## English

### Why

Hidradenitis suppurativa (HS, also called Verneuil's disease) is a chronic inflammatory skin disease. Adalimumab is the first approved biologic, but only part of the patients respond. Several public studies have measured gene expression in HS skin biopsies. This project asks a simple question: **can a model trained on some cohorts predict treatment response, or the molecular subtype, in a cohort it has never seen?**

### Data

Public GEO series of HS lesional skin, listed in [`configs/cohorts.yaml`](configs/cohorts.yaml). The list follows the reference study [*Classification of skin transcriptome reveals two molecular subtypes in hidradenitis suppurativa* (bioRxiv, 2025)](https://www.biorxiv.org/content/10.1101/2025.02.03.636243): 6 discovery cohorts (100 lesional samples) and 3 validation cohorts. Which cohorts carry treatment-response labels is still to be checked.

### Evaluation protocol

- **Leave-one-cohort-out:** every metric is computed on a cohort that was never used for training or preprocessing.
- **Per-cohort standardisation** first, to limit study effects without leaking information across cohorts.
- **Baseline:** elastic-net logistic regression, compared later with richer models.
- **Metric:** AUROC per held-out cohort, reported cohort by cohort, not only averaged.

### Layout

```
configs/cohorts.yaml     GEO cohorts and their role
src/hs_stratify/
  cohorts.py             cohort registry
  data.py                GEO download (GEOparse)
  preprocess.py          gene intersection and per-cohort standardisation
  splits.py              leave-one-cohort-out splits
  models.py              baseline model
  evaluate.py            cross-cohort evaluation
  cli.py                 command line
tests/                   unit tests on synthetic data
notebooks/               exploration
```

### Getting started

```bash
pip install -e ".[dev,data]"
pytest
hs-stratify cohorts          # list configured cohorts
hs-stratify download         # download all cohorts from GEO into data/raw
```

## Français

### Pourquoi

La maladie de Verneuil (hidradenitis suppurativa, HS) est une maladie inflammatoire chronique de la peau. L'adalimumab est le premier biomédicament autorisé, mais seule une partie des patients y répond. Plusieurs études publiques ont mesuré l'expression des gènes dans des biopsies de lésions. La question du projet : **un modèle entraîné sur certaines cohortes peut-il prédire la réponse au traitement, ou le sous-type moléculaire, dans une cohorte qu'il n'a jamais vue ?**

### Données

Des séries GEO publiques de peau lésionnelle, listées dans [`configs/cohorts.yaml`](configs/cohorts.yaml). La liste reprend l'étude de référence ([bioRxiv, 2025](https://www.biorxiv.org/content/10.1101/2025.02.03.636243)) : 6 cohortes de découverte (100 biopsies) et 3 cohortes de validation. Il reste à vérifier lesquelles contiennent la réponse au traitement.

### Protocole d'évaluation

- **Une cohorte laissée de côté à chaque fois :** toutes les métriques sont calculées sur une cohorte jamais utilisée pour l'entraînement ni le prétraitement.
- **Standardisation par cohorte** d'abord, pour limiter les effets d'étude sans faire fuiter d'information d'une cohorte à l'autre.
- **Baseline :** régression logistique elastic net, comparée ensuite à des modèles plus riches.
- **Métrique :** AUROC par cohorte testée, rapportée cohorte par cohorte et pas seulement en moyenne.

### Démarrer

```bash
pip install -e ".[dev,data]"
pytest
hs-stratify cohorts          # lister les cohortes configurées
hs-stratify download         # télécharger les cohortes depuis GEO dans data/raw
```

## License

MIT, see [LICENSE](LICENSE).
