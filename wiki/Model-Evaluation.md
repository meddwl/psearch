# Model Evaluation

PSearch performs two rounds of statistical evaluation: **internal** (during
model building, on training-set compounds) and **external** (after building,
on held-out test-set compounds). Both rounds use the same family of metrics.

---

## Confusion matrix terms

For each pharmacophore model screened against an evaluation set:

| Term | Meaning |
|------|---------|
| **TP** | Active compounds correctly identified as hits |
| **FP** | Inactive compounds incorrectly identified as hits |
| **P** | Total actives in the evaluation set |
| **N** | Total inactives in the evaluation set |
| **TN** | Inactives correctly not retrieved (N − FP) |

---

## Metrics

### Precision (positive predictive value)

```
Precision = TP / (TP + FP)
```

What fraction of the hits are truly active. High precision means few false
positives. This is the most important metric for prospective virtual screening,
where every selected compound is a candidate for synthesis or purchase.

### Recall (sensitivity, true positive rate)

```
Recall = TP / P
```

What fraction of all known actives the model retrieves. High recall means good
coverage of the active chemical space. Important when you want to find as many
active scaffolds as possible.

### F-scores

F-scores combine precision and recall into a single number. The parameter beta
controls the trade-off:

```
F_beta = (1 + beta^2) * (Precision * Recall) / (beta^2 * Precision + Recall)
```

| Score | beta | Emphasis | Used in PSearch for |
|-------|------|----------|---------------------|
| **F0.5** | 0.5 | Precision-weighted | Strategy 1 (centroid) model selection threshold: F0.5 >= 0.8 |
| **F1** | 1.0 | Balanced | Reported in `external_statistics.txt` |
| **F2** | 2.0 | Recall-weighted | Strategy 2 (per-cluster) model selection threshold: F2 >= 0.8 |

**Why F0.5 for Strategy 1?**
Centroid models are intended for prospective screening where false positives
waste experimental resources. F0.5 penalises false positives more heavily than
false negatives, so the selection threshold of F0.5 >= 0.8 keeps precision high
even if some actives are missed.

**Why F2 for Strategy 2?**
Per-cluster models aim to capture diverse binding modes across the active
chemical space. Missing an active (false negative) is a greater concern than
including a false positive, so recall is weighted more heavily. The F2 >= 0.8
threshold ensures broad coverage.

### Balanced accuracy (BA)

```
BA = (Recall + Specificity) / 2
     Specificity = TN / N
```

The average of the true positive rate and the true negative rate. BA = 0.5 is
equivalent to random guessing; BA = 1.0 is a perfect classifier. Reported in
`external_statistics.txt` but not used as a selection criterion.

### False positive rate (FPR)

```
FPR = FP / N
```

The fraction of inactives incorrectly flagged as hits. Reported alongside
recall so you can assess the recall/FPR trade-off in a ROC-style analysis.

### Enrichment factor (EF)

```
EF = Precision / (P / (P + N))
```

How many times richer in actives the hit list is compared to random selection
from the full compound set. EF = 1.0 means no enrichment; EF = 2.0 means the
hit list is twice as enriched in actives as a random sample. Useful for
comparing models when the active rate of your library is low.

In the example CDK8 dataset (94 actives, 133–134 inactives), top models reach
EF around 1.5, with precision around 0.5–0.6.

---

## Internal evaluation (during model building)

Computed inside `gen_pharm_models` on the training-set compounds themselves.
Used purely to prune the candidate model set at each feature-count step — a
model whose F0.5 (strategy 1) or F2 (strategy 2) falls below 0.8 is discarded
and not extended further.

Internal statistics are not a reliable estimate of prospective performance
because the same compounds are used for both model building and scoring. Use
them only to understand which feature counts produce viable models.

To save internal statistics for inspection:

```bash
gen_pharm_models -i dbs/cdk8.dat -ts trainset/t1.smi -o models/ -s
# writes models/intermediate_data/internal_statistics-t1-f4.txt, -f5.txt, ...
```

---

## External evaluation (after model building)

Computed by `external_stat` on the **external test set**: all compounds that
were **not** in the training set of a given model. Because each cluster's
training set contains only a subset of actives and inactives, the remainder
serve as an independent held-out test set. This is the reliable estimate of
prospective model performance.

```bash
external_stat -i cdk8.smi \
              -t my_project/trainset/ \
              -m my_project/models/ \
              -s my_project/raw_screen/ \
              -o my_project/external_statistics.txt
```

This step is run automatically as the final stage of `psearch`.

### Output file format

`external_statistics.txt` is a tab-separated file, one row per model, sorted
by recall descending then F0.5 descending:

```
target_id  model_id          TP  FP   P    N  precision  recall  FPR    F1     F2     F05    BA     EF     uniq_features  max_dist  features
cdk8       t2_f6_p1          83  79   94  133  0.512      0.883   0.594  0.648  0.771  0.559  0.644  1.237  3              5.121     aaAHH
cdk8       centroids_f5_p2   62  38   94  134  0.620      0.660   0.284  0.639  0.651  0.628  0.688  1.504  3              7.347     aDAH
```

The `features` column shows the feature composition; `uniq_features` counts
features with distinct 3D coordinates (a model with two features at the same
position counts as 1 unique feature). `max_dist` is the longest pairwise
distance between any two features in Å.

---

## Physical model descriptors

`get_pharm_props` computes structural properties for a set of model files,
useful for filtering or comparing models independently of activity data:

```bash
get_pharm_props -i my_project/models/ -o model_properties.txt
```

| Column | Meaning |
|--------|---------|
| `pharm_id` | Model filename without extension |
| `max_dist` | Maximum pairwise inter-feature distance (Å) |
| `nf` | Total feature count (including duplicates at the same position) |
| `nf_dist` | Features with unique 3D coordinates |
| `nf_polar` | Count of polar features (A, D, P, N) |
| `labels` | Concatenated feature-type string (e.g. `ADHH`) |

---

## Practical guidance

A model worth using prospectively typically has:

- **Precision >= 0.5** — at least half the hits are active
- **EF >= 1.2** — meaningfully better than random selection
- **Recall >= 0.3** — retrieves at least 30 % of known actives

The `features` column helps interpret the binding interaction captured:
`aDAH` describes a model with an aromatic centre, donor, acceptor, and
hydrophobic — suggestive of a mixed polar/hydrophobic binding pocket.
A model with only `aaaH` features describes a flat hydrophobic/aromatic
environment.
