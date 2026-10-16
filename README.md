# hs-stratify

Predicting adalimumab response in HS (Hidradenitis Suppurativa) from skin transcriptomics, with cross-cohort validation.

[English](#english) · [Français](#français)

> **Status:** work in progress. GEO data loading and the cross-cohort protocol are in place, with a first result (lesional vs non-lesional skin). Molecular subtypes and treatment response come next.

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

### First experiment: lesional vs non-lesional skin

Before predicting treatment response (only 26 labelled patients across two cohorts), the project checks that a model can tell **lesional** skin from **non-lesional or healthy** skin in a cohort it has never seen, including across technologies (RNA-seq to microarrays). Tissue labels come from each cohort's metadata (rules in `configs/cohorts.yaml`), using pre-treatment biopsies only; perilesional skin is left out.

Two normalisations are compared:
- **per cohort:** each gene centred and scaled within its cohort. This removes study effects, but also the lesion signal of a cohort made almost only of lesions;
- **per sample:** each sample centred and scaled over its own genes, independent of the other samples.

```bash
hs-stratify evaluate lesional   # AUROC per held-out cohort, written to results/lesional_loco.csv
```

**Results** (AUROC on each held-out cohort, 11,585 genes shared by the 9 cohorts):

| Held-out cohort | Technology | Samples (lesional) | Per sample | Per cohort |
|---|---|---|---|---|
| GSE154773 | RNA-seq | 32 (22) | 0.968 | 0.932 |
| GSE155176 | RNA-seq | 55 (26) | 1.000 | 1.000 |
| GSE189266 | RNA-seq | 18 (9) | 0.963 | 0.963 |
| GSE213761 | RNA-seq | 30 (15) | 0.862 | 0.853 |
| GSE249027 | RNA-seq | 27 (14) | 1.000 | 1.000 |
| GSE72702 | microarray | 30 (17) | 0.977 | 0.973 |
| GSE148027 | microarray | 22 (10) | 1.000 | 1.000 |
| GSE128637 | microarray | 21 (10) | 1.000 | 1.000 |
| **Mean** | | | **0.971** | **0.965** |

GSE151243 contains only lesions once perilesional skin is left out, so it is used for training but cannot be scored on its own.

What this shows:
- **The signal transfers across studies and technologies:** a model trained mostly on RNA-seq scores 0.97 to 1.00 on the three microarray cohorts it has never seen.
- **The two normalisations are close** (0.971 against 0.965). The gap is too small, with 8 cohorts, to claim one is better.
- **GSE213761 is the hardest cohort** (0.86). Its non-lesional biopsies come from patients with moderate to severe disease (Hurley II and III), whose "healthy-looking" skin may already carry inflammation.
- **This is a sanity check, not the end goal:** lesional skin is very different from non-lesional skin. The real question, predicting treatment response, is much harder and has far less data.

### Next: molecular subtypes

Lesional biopsies only (pre-treatment), centred within each cohort to remove study effects, on the 2,000 most variable genes:

1. **Consensus clustering** (Monti et al., 2003): k-means on 100 random subsets of 80% of the lesions, for k = 2 to 5. The PAC score (share of ambiguous pairs) picks k: close to 0 means clear-cut groups.
2. **Stability across cohorts:** for each cohort, subtypes are learnt without it, then its lesions are assigned to the nearest subtype. The adjusted Rand index (ARI) against the subtypes found on all data says whether they are reproducible (1 = identical, 0 = chance).
3. **Link with treatment response** on the 26 labelled patients (Fisher's exact test when k = 2), as an exploratory result.

On pure noise, this pipeline gives a PAC close to 1 and an ARI close to 0, so a high ARI on real data is not an artefact of the method.

```bash
hs-stratify evaluate subtypes   # tables written to results/subtypes_*.csv
```

### Layout

```
configs/cohorts.yaml     GEO cohorts and their role
src/hs_stratify/
  cohorts.py             cohort registry
  geo.py                 GEO URLs and file parsers
  data.py                builds one cohort's expression matrix
  labels.py              tissue and treatment response labels (field or HiSCR)
  subtypes.py            consensus clustering and cross-cohort stability
  experiments.py         reproducible experiments
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

### Première expérience : peau lésionnelle ou non

Avant de prédire la réponse au traitement (26 patients étiquetés seulement, sur deux cohortes), le projet vérifie qu'un modèle sait distinguer une peau **lésionnelle** d'une peau **non lésionnelle ou saine** dans une cohorte qu'il n'a jamais vue, y compris en changeant de technologie (du RNA-seq aux puces). Le type de tissu vient des métadonnées de chaque cohorte (règles dans `configs/cohorts.yaml`), sur les biopsies avant traitement uniquement ; la peau péri-lésionnelle est écartée.

Deux normalisations sont comparées :
- **par cohorte :** chaque gène est centré-réduit au sein de sa cohorte. Cela retire les effets d'étude, mais aussi le signal de lésion d'une cohorte composée presque uniquement de lésions ;
- **par échantillon :** chaque échantillon est centré-réduit sur ses propres gènes, sans dépendre des autres échantillons.

```bash
hs-stratify evaluate lesional   # AUROC par cohorte laissée de côté, écrite dans results/lesional_loco.csv
```

**Résultats** (AUROC sur chaque cohorte laissée de côté, 11 585 gènes communs aux 9 cohortes) :

| Cohorte testée | Technologie | Échantillons (lésions) | Par échantillon | Par cohorte |
|---|---|---|---|---|
| GSE154773 | RNA-seq | 32 (22) | 0,968 | 0,932 |
| GSE155176 | RNA-seq | 55 (26) | 1,000 | 1,000 |
| GSE189266 | RNA-seq | 18 (9) | 0,963 | 0,963 |
| GSE213761 | RNA-seq | 30 (15) | 0,862 | 0,853 |
| GSE249027 | RNA-seq | 27 (14) | 1,000 | 1,000 |
| GSE72702 | puce | 30 (17) | 0,977 | 0,973 |
| GSE148027 | puce | 22 (10) | 1,000 | 1,000 |
| GSE128637 | puce | 21 (10) | 1,000 | 1,000 |
| **Moyenne** | | | **0,971** | **0,965** |

GSE151243 ne contient que des lésions une fois la peau péri-lésionnelle écartée : elle sert à l'entraînement mais ne peut pas être évaluée seule.

Ce que ça montre :
- **Le signal se transfère d'une étude et d'une technologie à l'autre :** un modèle entraîné surtout sur du RNA-seq obtient de 0,97 à 1,00 sur les trois cohortes de puces qu'il n'a jamais vues.
- **Les deux normalisations sont proches** (0,971 contre 0,965). Avec 8 cohortes, l'écart est trop faible pour affirmer que l'une est meilleure.
- **GSE213761 est la cohorte la plus difficile** (0,86). Ses biopsies non lésionnelles viennent de patients à un stade modéré ou sévère (Hurley II et III), dont la peau d'apparence saine porte peut-être déjà de l'inflammation.
- **C'est un contrôle de cohérence, pas l'objectif final :** une peau lésionnelle est très différente d'une peau non lésionnelle. La vraie question, prédire la réponse au traitement, est bien plus difficile et dispose de beaucoup moins de données.

### Suite : sous-types moléculaires

Biopsies de lésion uniquement (avant traitement), centrées au sein de chaque cohorte pour retirer les effets d'étude, sur les 2 000 gènes les plus variables :

1. **Clustering de consensus** (Monti et al., 2003) : k-means sur 100 sous-ensembles aléatoires de 80 % des lésions, pour k = 2 à 5. Le score PAC (part des paires ambiguës) choisit k : proche de 0, les groupes sont nets.
2. **Stabilité entre cohortes :** pour chaque cohorte, les sous-types sont appris sans elle, puis ses lésions sont rattachées au sous-type le plus proche. L'indice de Rand ajusté (ARI) par rapport aux sous-types trouvés sur toutes les données dit s'ils sont reproductibles (1 = identiques, 0 = hasard).
3. **Lien avec la réponse au traitement** sur les 26 patients étiquetés (test exact de Fisher si k = 2), à titre exploratoire.

Sur du bruit pur, ce pipeline donne un PAC proche de 1 et un ARI proche de 0 : un ARI élevé sur les vraies données n'est donc pas un artefact de la méthode.

```bash
hs-stratify evaluate subtypes   # tableaux écrits dans results/subtypes_*.csv
```

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
