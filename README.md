# hs-stratify

Predicting adalimumab response in HS (Hidradenitis Suppurativa) from skin transcriptomics, with cross-cohort validation.

[English](#english) · [Français](#français)

> **Status:** complete. Three experiments on 9 public cohorts (435 samples): lesional skin is recognised across studies, two molecular subtypes are reproduced, and treatment response cannot be predicted from the 26 labelled patients available. See [Conclusions](#conclusions).

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

### Second experiment: molecular subtypes

Lesional biopsies only (pre-treatment), centred within each cohort to remove study effects, on the 2,000 most variable genes:

1. **Consensus clustering** (Monti et al., 2003): k-means on 100 random subsets of 80% of the lesions, for k = 2 to 5. The PAC score (share of ambiguous pairs) picks k: close to 0 means clear-cut groups.
2. **Stability across cohorts:** for each cohort, subtypes are learnt without it, then its lesions are assigned to the nearest subtype. The adjusted Rand index (ARI) against the subtypes found on all data says whether they are reproducible (1 = identical, 0 = chance).
3. **Link with treatment response** on the 26 labelled patients (Fisher's exact test when k = 2), as an exploratory result.

On pure noise, this pipeline gives a PAC close to 1 and an ARI close to 0, so a high ARI on real data is not an artefact of the method.

```bash
hs-stratify evaluate subtypes   # tables written to results/subtypes_*.csv
```

**Results** (143 pre-treatment lesions, 11,585 shared genes):

- **Two subtypes** (PAC 0.23 for k = 2, against 0.47 to 0.49 for k = 3 to 5). The groups are reasonably clear but not clear-cut, which fits a gradient between two biological programmes more than two separate diseases.
- **They are reproducible across cohorts:** learnt without a cohort, they are found again on it with a pooled ARI of 0.84 (96% agreement). Per cohort, ARI ranges from 0.56 to 1.00, microarray cohorts included. Both subtypes appear in every cohort, but part of this balance comes from centring each cohort.
- **They match the reference study's two subtypes:**
  - subtype 0, **immune / inflammatory**: leukocyte genes (SELL, PLAC8, IL10RA, MS4A7, HCLS1, LY96, WAS, SASH3);
  - subtype 1, **epidermal / keratinisation**: keratinocyte differentiation genes (DLX3, PKP1, KCNK7, GGT6, POU3F1, BMP7, AHNAK2).
- **No link with treatment response here:** 60% responders in subtype 0 against 64% in subtype 1 (26 patients, Fisher p = 1.0). The reference study reports 62.5% against 20%. With 26 patients this test can only detect a very large effect, and response labels were defined differently (anti-TNF response field for GSE155176, HiSCR recomputed for GSE213761), so this does not refute it, but the association is not reproduced.

### Third experiment: predicting treatment response

Only two cohorts carry response labels: GSE155176 (9 responders, 3 non-responders) and GSE213761 (7 and 7). With so few patients a single AUROC means little, so each score comes with:
- a **95% bootstrap interval**, resampling responders and non-responders separately;
- a **permutation test**: the model is retrained on shuffled training labels, and the p-value is the share of shuffled runs that do as well as the real one.

Settings:
- **cross-cohort**, the honest one: train on one cohort, test on the other, in both directions;
- **pooled leave-one-out**, for comparison: each patient is scored by a model trained on the 25 others, both cohorts mixed. On pure noise, leave-one-out AUROC sits well below 0.5 (a known bias with tiny samples), which is why the permutation p-value, not the AUROC, is the number to read.

The model is an L2 logistic regression on the 2,000 most variable genes of the training set: the permutation test refits it thousands of times, which the elastic net is too slow for on a laptop. On pure noise, the cross-cohort test gave no p-value under 0.05 in 20 simulated runs.

```bash
hs-stratify evaluate response   # written to results/response.csv
```

**Results** (26 patients, 17,026 genes shared by the two cohorts):

| Setting | Train | Test | AUROC | 95% interval | Permutation p |
|---|---|---|---|---|---|
| Cross-cohort | GSE155176 (9R/3NR) | GSE213761 (7R/7NR) | 0.49 | 0.16 to 0.80 | 0.54 |
| Cross-cohort | GSE213761 (7R/7NR) | GSE155176 (9R/3NR) | 0.70 | 0.44 to 1.00 | 0.21 |
| Pooled leave-one-out | 25 patients | each patient | 0.33 | 0.13 to 0.53 | 0.78 |

**No setting beats chance.** The 0.70 in one direction is tempting, but its interval goes down to 0.44 and shuffled labels do as well one time in five. The pooled leave-one-out at 0.33 is the small-sample bias described above, not a model that gets things backwards.

This is a sample size problem before being a model problem. By Hanley and McNeil's formula, estimating an AUROC of 0.75 to within ±0.1 takes about 40 responders and 40 non-responders, three times the patients available here, with the same response definition across studies.

## Conclusions

1. **Public HS transcriptomes can be pooled:** after per-sample normalisation, a model trained on some studies recognises lesional skin in a study it has never seen (mean AUROC 0.97), including from RNA-seq to microarrays.
2. **The two molecular subtypes of the reference study are reproduced** (immune vs epidermal) and hold up across cohorts (ARI 0.84), with an independent pipeline.
3. **Predicting adalimumab response is not possible with public data today:** neither the subtypes nor a transcriptome-wide model do better than chance on 26 patients, and the subtype–response link reported by the reference study is not reproduced. Answering that question needs a cohort of about 80 treated patients.

### Layout

```
configs/cohorts.yaml     GEO cohorts and their role
src/hs_stratify/
  cohorts.py             cohort registry
  geo.py                 GEO URLs and file parsers
  data.py                builds one cohort's expression matrix
  labels.py              tissue and treatment response labels (field or HiSCR)
  subtypes.py            consensus clustering and cross-cohort stability
  response.py            bootstrap intervals and permutation tests for tiny samples
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

### Deuxième expérience : sous-types moléculaires

Biopsies de lésion uniquement (avant traitement), centrées au sein de chaque cohorte pour retirer les effets d'étude, sur les 2 000 gènes les plus variables :

1. **Clustering de consensus** (Monti et al., 2003) : k-means sur 100 sous-ensembles aléatoires de 80 % des lésions, pour k = 2 à 5. Le score PAC (part des paires ambiguës) choisit k : proche de 0, les groupes sont nets.
2. **Stabilité entre cohortes :** pour chaque cohorte, les sous-types sont appris sans elle, puis ses lésions sont rattachées au sous-type le plus proche. L'indice de Rand ajusté (ARI) par rapport aux sous-types trouvés sur toutes les données dit s'ils sont reproductibles (1 = identiques, 0 = hasard).
3. **Lien avec la réponse au traitement** sur les 26 patients étiquetés (test exact de Fisher si k = 2), à titre exploratoire.

Sur du bruit pur, ce pipeline donne un PAC proche de 1 et un ARI proche de 0 : un ARI élevé sur les vraies données n'est donc pas un artefact de la méthode.

```bash
hs-stratify evaluate subtypes   # tableaux écrits dans results/subtypes_*.csv
```

**Résultats** (143 lésions avant traitement, 11 585 gènes communs) :

- **Deux sous-types** (PAC de 0,23 pour k = 2, contre 0,47 à 0,49 pour k = 3 à 5). Les groupes sont assez nets sans être tranchés, ce qui évoque un gradient entre deux programmes biologiques plutôt que deux maladies distinctes.
- **Ils sont reproductibles d'une cohorte à l'autre :** appris sans une cohorte, ils s'y retrouvent avec un ARI global de 0,84 (96 % d'accord). Par cohorte, l'ARI va de 0,56 à 1,00, cohortes de puces comprises. Les deux sous-types sont présents dans toutes les cohortes, mais ce bon équilibre vient en partie du centrage par cohorte.
- **Ils correspondent aux deux sous-types de l'étude de référence :**
  - sous-type 0, **immunitaire / inflammatoire** : gènes des globules blancs (SELL, PLAC8, IL10RA, MS4A7, HCLS1, LY96, WAS, SASH3) ;
  - sous-type 1, **épidermique / kératinisation** : gènes de différenciation des kératinocytes (DLX3, PKP1, KCNK7, GGT6, POU3F1, BMP7, AHNAK2).
- **Pas de lien avec la réponse au traitement ici :** 60 % de répondeurs dans le sous-type 0 contre 64 % dans le sous-type 1 (26 patients, Fisher p = 1,0). L'étude de référence rapporte 62,5 % contre 20 %. Avec 26 patients, ce test ne peut détecter qu'un effet très fort, et la réponse n'est pas définie de la même façon (champ de réponse aux anti-TNF pour GSE155176, HiSCR recalculé pour GSE213761) : cela ne la contredit pas, mais l'association n'est pas reproduite.

### Troisième expérience : prédire la réponse au traitement

Seules deux cohortes ont la réponse au traitement : GSE155176 (9 répondeurs, 3 non-répondeurs) et GSE213761 (7 et 7). Avec si peu de patients, une AUROC seule ne veut pas dire grand-chose : chaque score est donc accompagné
- d'un **intervalle à 95 % par bootstrap**, en rééchantillonnant répondeurs et non-répondeurs séparément ;
- d'un **test de permutation** : le modèle est réentraîné sur des étiquettes mélangées, et la p-value est la part de ces essais qui font aussi bien que le vrai.

Deux cadres :
- **d'une cohorte à l'autre**, le cadre honnête : entraînement sur une cohorte, test sur l'autre, dans les deux sens ;
- **leave-one-out groupé**, pour comparaison : chaque patient est évalué par un modèle entraîné sur les 25 autres, les deux cohortes mélangées. Sur du bruit pur, l'AUROC en leave-one-out tombe nettement sous 0,5 (un biais connu sur les tout petits échantillons) : c'est donc la p-value de permutation, pas l'AUROC, qu'il faut lire.

Le modèle est une régression logistique L2 sur les 2 000 gènes les plus variables de l'entraînement : le test de permutation le réentraîne des milliers de fois, ce qui serait trop lent avec l'elastic net sur un portable. Sur du bruit pur, le test d'une cohorte à l'autre n'a donné aucune p-value sous 0,05 en 20 simulations.

```bash
hs-stratify evaluate response   # écrit dans results/response.csv
```

**Résultats** (26 patients, 17 026 gènes communs aux deux cohortes) :

| Cadre | Entraînement | Test | AUROC | Intervalle à 95 % | p de permutation |
|---|---|---|---|---|---|
| D'une cohorte à l'autre | GSE155176 (9R/3NR) | GSE213761 (7R/7NR) | 0,49 | 0,16 à 0,80 | 0,54 |
| D'une cohorte à l'autre | GSE213761 (7R/7NR) | GSE155176 (9R/3NR) | 0,70 | 0,44 à 1,00 | 0,21 |
| Leave-one-out groupé | 25 patients | chaque patient | 0,33 | 0,13 à 0,53 | 0,78 |

**Aucun cadre ne fait mieux que le hasard.** Le 0,70 dans un sens est tentant, mais son intervalle descend à 0,44 et des étiquettes mélangées font aussi bien une fois sur cinq. Le 0,33 du leave-one-out groupé est le biais des petits échantillons décrit plus haut, pas un modèle qui se tromperait systématiquement.

C'est un problème d'effectif avant d'être un problème de modèle. D'après la formule de Hanley et McNeil, estimer une AUROC de 0,75 à ±0,1 près demande environ 40 répondeurs et 40 non-répondeurs, soit trois fois plus de patients qu'ici, avec la même définition de la réponse d'une étude à l'autre.

## Conclusions

1. **Les transcriptomes publics de la maladie de Verneuil peuvent être réunis :** après normalisation par échantillon, un modèle entraîné sur certaines études reconnaît la peau lésionnelle dans une étude qu'il n'a jamais vue (AUROC moyenne de 0,97), y compris en passant du RNA-seq aux puces.
2. **Les deux sous-types moléculaires de l'étude de référence sont retrouvés** (immunitaire et épidermique) et tiennent d'une cohorte à l'autre (ARI de 0,84), avec un pipeline indépendant.
3. **Prédire la réponse à l'adalimumab n'est pas possible aujourd'hui avec les données publiques :** ni les sous-types ni un modèle sur tout le transcriptome ne font mieux que le hasard sur 26 patients, et le lien entre sous-type et réponse rapporté par l'étude de référence n'est pas reproduit. Répondre à cette question demanderait une cohorte d'environ 80 patients traités.

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
