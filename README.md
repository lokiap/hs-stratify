# hs-stratify

Predicting adalimumab response in HS (Hidradenitis Suppurativa) from skin transcriptomics, with cross-cohort validation.

[English](#english) · [Français](#français)

> **Status:** work in progress. Cohort registry, GEO data loading and the evaluation protocol are in place; models and results are not there yet.

## English

### Why

Hidradenitis suppurativa (HS, also called Verneuil's disease) is a chronic inflammatory skin disease. Adalimumab is the first approved biologic, but only part of the patients respond. Several public studies have measured gene expression in HS skin biopsies. This project asks a simple question: **can a model trained on some cohorts predict treatment response, or the molecular subtype, in a cohort it has never seen?**

### Data

Public GEO series of HS lesional skin, listed in [`configs/cohorts.yaml`](configs/cohorts.yaml). The list follows the reference study [*Classification of skin transcriptome reveals two molecular subtypes in hidradenitis suppurativa* (bioRxiv, 2025)](https://www.biorxiv.org/content/10.1101/2025.02.03.636243): 6 discovery cohorts (RNA-seq, 100 lesional samples) and 3 validation cohorts (microarrays). Two cohorts carry treatment response: GSE155176 gives it directly (anti-TNF, baseline lesional biopsies), and for GSE213761 (adalimumab) it is recomputed as HiSCR from abscess, nodule and fistula counts between visits V1 and V12.

- **RNA-seq:** raw counts recomputed by NCBI on GRCh38, the same pipeline for every study, fetched from the GEO FTP (the website's download page can answer a script with a reCAPTCHA), then log2(CPM + 1) with the reference study's filter (CPM ≥ 1 in ≥ 20% of samples).
- **Microarrays:** values from the GEO series matrix, probes mapped to genes with the GPL annotation and averaged per gene.
- **Sample metadata** (tissue, response, timepoint) comes from the series matrix of each cohort.

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
  geo.py                 GEO URLs and file parsers
  data.py                builds one cohort's expression matrix
  labels.py              treatment response labels (field or HiSCR)
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
pip install -e ".[dev]"
pytest
hs-stratify cohorts          # list configured cohorts
hs-stratify inspect GSE155176  # show a cohort's metadata fields
hs-stratify build            # download from GEO and write data/processed/
```

## Français

### Pourquoi

La maladie de Verneuil (hidradenitis suppurativa, HS) est une maladie inflammatoire chronique de la peau. L'adalimumab est le premier biomédicament autorisé, mais seule une partie des patients y répond. Plusieurs études publiques ont mesuré l'expression des gènes dans des biopsies de lésions. La question du projet : **un modèle entraîné sur certaines cohortes peut-il prédire la réponse au traitement, ou le sous-type moléculaire, dans une cohorte qu'il n'a jamais vue ?**

### Données

Des séries GEO publiques de peau lésionnelle, listées dans [`configs/cohorts.yaml`](configs/cohorts.yaml). La liste reprend l'étude de référence ([bioRxiv, 2025](https://www.biorxiv.org/content/10.1101/2025.02.03.636243)) : 6 cohortes de découverte (RNA-seq, 100 biopsies) et 3 cohortes de validation (puces). Deux cohortes contiennent la réponse au traitement : GSE155176 la donne directement (anti-TNF, biopsies de lésion avant traitement), et pour GSE213761 (adalimumab) elle est recalculée selon le critère HiSCR à partir des comptes d'abcès, de nodules et de fistules entre les visites V1 et V12.

- **RNA-seq :** comptages recalculés par NCBI sur GRCh38, le même traitement pour toutes les études, récupérés sur le FTP de GEO (la page de téléchargement du site web peut renvoyer un reCAPTCHA à un script), puis log2(CPM + 1) avec le filtre de l'étude de référence (CPM ≥ 1 dans au moins 20 % des échantillons).
- **Puces :** valeurs du series matrix GEO, sondes rattachées aux gènes via l'annotation GPL et moyennées par gène.
- **Métadonnées des échantillons** (tissu, réponse, temps) : lues dans le series matrix de chaque cohorte.

### Protocole d'évaluation

- **Une cohorte laissée de côté à chaque fois :** toutes les métriques sont calculées sur une cohorte jamais utilisée pour l'entraînement ni le prétraitement.
- **Standardisation par cohorte** d'abord, pour limiter les effets d'étude sans faire fuiter d'information d'une cohorte à l'autre.
- **Baseline :** régression logistique elastic net, comparée ensuite à des modèles plus riches.
- **Métrique :** AUROC par cohorte testée, rapportée cohorte par cohorte et pas seulement en moyenne.

### Démarrer

```bash
pip install -e ".[dev]"
pytest
hs-stratify cohorts          # lister les cohortes configurées
hs-stratify inspect GSE155176  # afficher les champs de métadonnées d'une cohorte
hs-stratify build            # télécharger depuis GEO et écrire data/processed/
```

## License

MIT, see [LICENSE](LICENSE).
