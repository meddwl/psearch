# Activity Prediction

## Overview

After screening a compound library against pharmacophore models, the
`prediction` tool converts hit/miss results into a **probability of activity**
for each compound towards each target. It does this by weighting each model's
hit indicator by that model's empirical **precision** — the fraction of its
historical hits that were truly active — and then combining scores across all
models for the same target via a consensus scheme.

This approach is described in full in:

> Probabilistic Approach for Virtual Screening Based on Multiple Pharmacophores
> Madzhidov et al., *Molecules* **2020**, 25(2), 385
> https://doi.org/10.3390/molecules25020385

---

## Inputs required

### 1. VS results directory

The output directory of `screen_db`, containing one `.txt` hit-list file per
model. File names must follow the convention `<target_id>.<model_id>.txt`
(e.g. `cdk8.t1_f7_p3.txt`). This naming is produced automatically by
`screen_db`.

### 2. Model precision file

A tab-separated file with three required columns: `target_id`, `model_id`,
`precision`.

- **After running your own models**: pass `external_statistics.txt` (produced
  by `external_stat` or the `psearch` pipeline). The `precision` column there
  is computed on the held-out external test set, giving an unbiased estimate.
- **Multi-target profiling with bundled ChEMBL models**: omit `-p` entirely.
  PSearch uses a pre-computed precision file bundled at
  `psearch/pharmacophores/pharmacophores_stat.csv`.

---

## Running prediction

### After your own model-building pipeline

```bash
prediction -s my_project/vs/models/ \
           -p my_project/external_statistics.txt \
           -f max \
           -o my_project/results.txt
```

### Multi-target profiling with bundled ChEMBL models

```bash
prediction -s profiling/vs/ -o profiling/result_multiprofiling.txt
```

Omitting `-p` uses the pre-computed precision statistics for the bundled
ChEMBL pharmacophore library.

---

## Scoring schemes

The `-f/--scoring_scheme` flag controls how individual model probabilities are
combined into a single per-target score:

| Scheme | Behaviour | When to use |
|--------|-----------|-------------|
| `max` | Score = precision of the highest-precision model that hit the compound | Compound is predicted active if **any** high-precision model matches. Conservative; recommended for prospective screening. |
| `mean` | Score = average precision across all models that hit the compound | Aggregates evidence from multiple models. More informative when many models contribute, as in multi-target profiling. |

If a compound is not hit by any model for a target, its score for that target
is `NaN` — no prediction is made.

---

## Output format

The output is a tab-separated file. Rows are compounds; columns are targets.
Values are activity probability scores rounded to three decimal places. Rows
are sorted by the first target column in descending order.

```
mol_id          cdk8    CHEMBL1978    CHEMBL205
CHEMBL3798663   0.883   NaN           0.512
CHEMBL3800311   0.620   0.495         NaN
CHEMBL3798944   0.512   NaN           NaN
```

`NaN` means no pharmacophore model for that target matched the compound.

---

## Interpreting scores

The score represents the empirical precision of the model (or models) that
matched the compound:

- A score of **0.88** means the model(s) that matched this compound
  historically had 88 % of their hits confirmed as actives in the external
  validation set.
- A score of **0.50** means roughly half of hits from matching models were
  active. Still a useful enrichment when the baseline active rate in a random
  library sample is much lower than 50 %, but weaker evidence than 0.88.

There is no single universal cutoff. The appropriate threshold depends on:
- The baseline active rate of your compound library
- How many false positives you can tolerate
- Whether you are doing prospective screening (prefer high cutoff) or
  hypothesis generation (can tolerate lower cutoff)

The practical approach is to sort by score, inspect the top-ranked compounds,
and set a cutoff based on how many compounds you can follow up experimentally.

---

## Using precision from your own models

Run `external_stat` after model building to compute precision on the held-out
test set, then pass the result to `prediction`:

```bash
# Step 1: compute external statistics (run automatically by psearch pipeline)
external_stat -i cdk8.smi \
              -t my_project/trainset/ \
              -m my_project/models/ \
              -s my_project/raw_screen/ \
              -o my_project/external_statistics.txt

# Step 2: screen a new library
screen_db -d dbs/new_library.dat \
          -q my_project/models/ \
          -o my_project/vs/ -c 4

# Step 3: predict activity probabilities
prediction -s my_project/vs/ \
           -p my_project/external_statistics.txt \
           -f max \
           -o my_project/results.txt
```

The `precision` column in `external_statistics.txt` is computed on the external
test set (compounds not in the training set of each model), so it is an
unbiased estimate of how the model will perform on unseen compounds.

---

## Full multi-target profiling workflow

This workflow screens a compound library against all bundled ChEMBL models and
ranks compounds by predicted activity across all targets:

```bash
# Build the database for your molecules of interest
gen_db -i mols_for_profiling.smi -d dbs/mols_for_profiling.dat -c 4 -v

# Screen against all bundled ChEMBL pharmacophores (no -q flag)
screen_db -d dbs/mols_for_profiling.dat -o multiprofiling/vs/ -c 4 -v

# Compute activity probabilities across all targets
prediction -s multiprofiling/vs/ -o multiprofiling/result_multiprofiling.txt
```

The output lists each compound's predicted probability of activity for each
ChEMBL target represented in the bundled model library.
