# Ligand-Based 3D Pharmacophore Modeling

## Overview

PSearch constructs pharmacophore models from a labelled compound set — molecules
annotated as active (1) or inactive (0) against a target. The process has three
stages:

1. **Database generation** — enumerate stereoisomers, embed 3D conformers, and
   compute pharmacophore feature coordinates for every compound.
2. **Training set selection** — cluster compounds by chemical similarity and
   select representative subsets as training sets.
3. **Iterative model building** — enumerate 4-feature pharmacophore sub-graphs
   from training-set conformers, score them by internal precision/recall
   statistics, then repeatedly extend surviving models by one feature at a time
   until quality can no longer improve.

## Pharmacophore feature types

PSearch uses the pmapper library's default feature definitions. Six feature
types are recognised:

| Code | Meaning |
|------|---------|
| `a` | Aromatic ring centre |
| `A` | Hydrogen-bond acceptor |
| `D` | Hydrogen-bond donor |
| `H` | Hydrophobic centre |
| `P` | Positively charged group |
| `N` | Negatively charged group |

A pharmacophore model is a set of these features with defined 3D coordinates.
The `features` column of `external_statistics.txt` shows each model's feature
composition (e.g. `aaAHH` = two aromatic centres + one acceptor + two
hydrophobics).

---

## Step 1 — Build the database

The database pre-computes stereoisomers, conformers, and pharmacophore
fingerprints so that model building and screening do not need to repeat
expensive 3D generation.

```bash
gen_db -i cdk8.smi -d dbs/cdk8.dat -c 4 -v
```

### Input file format

A tab-separated SMILES file with (at minimum) two columns:

```
SMILES<TAB>mol_id<TAB>activity
O=C1NCCC12CCN(...)CC2    CHEMBL3798663    1
c1ccc2ccccc2c1           CHEMBL000001     0
```

The activity column (1/0) is required for model building but not for
screening-only workflows.

### Key parameters

| Flag | Default | What it controls |
|------|---------|-----------------|
| `-n/--nconf` | 50 | Conformers per stereoisomer. More conformers improve pharmacophore coverage but increase database size and build time. |
| `-s/--nstereo` | 5 | Max stereoisomers per compound. Stereocentres with explicit configuration are preserved. |
| `-e/--energy_cutoff` | None | Discard conformers with MMFF energy > N kcal/mol above the lowest-energy conformer. |
| `-r/--rms` | None | Discard conformers with pairwise RMSD < N Å (removes near-duplicate geometries). |
| `-b/--bin_step` | 1 | Coordinate binning resolution in Å. **Must be the same value for all subsequent steps.** |
| `--seed` | -1 | Random seed for reproducible conformer generation. -1 means no fixed seed. |

### Database files

The database is stored as two sidecar files: `cdk8.dat` and `cdk8.dir`. Both
are required; never move one without the other.

### Duplicate handling

- If the same SMILES appears twice, the second occurrence is skipped with a
  warning.
- If the same `mol_id` appears for two distinct structures, a `#N` suffix is
  appended to the later entry.
- When any duplicates or renames occur, a corrected `*-updated.smi` file is
  written alongside the original input.

---

## Step 2 — Training set selection

Compounds are clustered by Tanimoto similarity using the Butina algorithm.
Two fingerprint types are available:

- **2D pharmacophore fingerprints** (default) — based on feature-pair distances
- **FCFP4** (add `--fcfp4`) — Morgan radius-2 feature-based fingerprints

The Butina distance cutoff (`-t`, default 0.4) controls cluster granularity:
lower values create more, smaller clusters and therefore more per-cluster
training sets.

Training set selection can be run as part of the full pipeline (`psearch`) or
as a standalone step:

```bash
select_training_set -i cdk8.smi -o my_project/trainset/ -ts 1 2 -t 0.4
```

### Strategy 1 — Centroid training set

One training set is formed from the **centroids** (most representative members)
of every cluster that contains at least 5 compounds. Actives and inactives are
clustered separately so the centroid set spans chemical diversity. This produces
a **single file**: `trainset/centroids.smi`.

**Best for**: prospective virtual screening where false positives are costly.
Models built from centroid training sets are precision-biased (selected at
F0.5 >= 0.8 during internal evaluation).

### Strategy 2 — Per-cluster training sets

One training set is created **per cluster** that contains at least 5 actives.
Each training set contains up to 5 actives from that cluster, up to 5 inactives
from that cluster, and the inactive centroids from all clusters (to maintain
broad inactive coverage). This produces multiple files: `t0.smi`, `t1.smi`, …
`tN.smi`.

**Best for**: retrospective validation studies or when you suspect the target
has multiple binding modes (e.g. allosteric vs. orthosteric). Models are
recall-biased (selected at F2 >= 0.8 during internal evaluation), aiming for
broad active-space coverage.

### Running both strategies

Both strategies run by default:

```bash
psearch -i cdk8.smi -d dbs/cdk8.dat -p my_project/ -m 1 2  # default
psearch -i cdk8.smi -d dbs/cdk8.dat -p my_project/ -m 1     # centroid only
psearch -i cdk8.smi -d dbs/cdk8.dat -p my_project/ -m 2     # per-cluster only
```

---

## Step 3 — Iterative model building

For each training set, `gen_pharm_models` carries out the following loop:

1. Enumerate all **4-feature sub-graphs** from the pharmacophore of every
   active conformer in the training set.
2. Compute internal TP, FP, precision, recall, F2, and F0.5 for each unique
   pharmacophore hash.
3. Prune to the top-scoring models (F0.5 >= 0.8 for strategy 1; F2 >= 0.8 for
   strategy 2).
4. **Extend** each surviving model by one additional feature from its source
   conformer.
5. Repeat from step 2 until the feature limit (`--upper`) is reached or adding
   features no longer yields any models above threshold.

This can be run as part of `psearch` or standalone:

```bash
gen_pharm_models -i dbs/cdk8.dat -ts my_project/trainset/t1.smi \
    -o my_project/models/ -b 1
```

The `--bin_step` value **must match** the value used in `gen_db`.

### Controlling model complexity

| Flag | Default | Effect |
|------|---------|--------|
| `-l/--lower` | 3 | Starting feature count for model generation |
| `-u/--upper` | None | Maximum feature count; None = no limit |
| `-f/--save_model_complexity` | None | Save models at all feature counts >= this value; None = save only the final set |

By default, only the most complex (highest-feature-count) surviving models are
saved. To also retain intermediate-complexity models:

```bash
gen_pharm_models -i dbs/cdk8.dat -ts trainset/t1.smi -o models/ -f 4
# saves models at 4 features, 5 features, ..., and final complexity
```

### Inspecting intermediate statistics

Add `-s/--save_statistics` to write per-step internal statistics files to
`models/intermediate_data/`:

```bash
gen_pharm_models -i dbs/cdk8.dat -ts trainset/t1.smi -o models/ -s
```

---

## Model file format and naming convention

Each pharmacophore model is stored as a plain `.xyz` text file. The filename
encodes its provenance:

```
<database_name>.<cluster_id>_f<num_features>_p<index>.xyz
```

Examples:

```
cdk8.centroids_f5_p0.xyz   — centroid strategy, 5 features, model index 0
cdk8.t1_f7_p3.xyz          — cluster 1 training set, 7 features, model index 3
```

The file content is feature labels and 3D coordinates in XYZ format:

```
A   1.234   5.678  -0.321
H   3.100   2.200   1.500
D  -0.500   4.100   2.300
```

Models are portable plain-text files and can be shared, inspected, or loaded
without access to the original database:

```python
from pmapper.pharmacophore import Pharmacophore
p = Pharmacophore()
p.load_from_xyz("cdk8.t1_f7_p3.xyz")
print(p.get_feature_coords())
```

---

## Bundled ChEMBL models

PSearch ships with pre-built pharmacophore models for numerous ChEMBL targets
in `psearch/pharmacophores/chembl_models/`. These follow the same naming
convention (`CHEMBL<id>.<cluster_id>_pharm<N>_<index>.xyz`) and are used
automatically by `screen_db` when no `-q` argument is supplied.
