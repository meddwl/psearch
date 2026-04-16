# PSearch Wiki

## What PSearch does

PSearch generates 3D ligand-based pharmacophore models from a set of active
and inactive compounds and uses those models to screen compound databases. A
pharmacophore describes the spatial arrangement of chemical features —
hydrogen-bond donors and acceptors, hydrophobic centres, charged groups, and
aromatic rings — that are necessary for a molecule to bind a target.

The key idea is that a pharmacophore is derived computationally from 3D
conformers of known active compounds, filtered so that the resulting model
is selective (few false positives) while retaining broad coverage (few false
negatives). PSearch iteratively builds candidate models starting from 4-feature
sub-graphs and expands them one feature at a time, keeping only those that pass
internal precision/recall thresholds.

PSearch can also run in **multi-target profiling mode** using a bundled library
of models built from ChEMBL data, without requiring any user-supplied training
compounds.

## When to use PSearch

| Use case | Recommended workflow |
|----------|---------------------|
| You have actives and inactives for one target and want hit-finding models | Full pipeline: `gen_db` → `psearch` → `screen_db` |
| You want to screen a compound database against your own models | `gen_db` → `screen_db` |
| You want to profile molecules against many ChEMBL targets at once | `gen_db` → `screen_db` (no `-q` flag) → `prediction` |
| You want to validate models rigorously before prospective use | Full pipeline includes automatic external validation |

PSearch suits both retrospective validation studies and prospective virtual
screening campaigns. Strategy 1 (centroid models) produces high-precision
models suitable for prospective screening. Strategy 2 (per-cluster models)
gives better recall, which is more appropriate when assessing coverage across
diverse binding modes.

## Quick start

```bash
# 1. Build the conformer/pharmacophore database from your labelled SMILES file
gen_db -i cdk8.smi -d dbs/cdk8.dat -c 4 -v

# 2. Run the full model-building and validation pipeline
psearch -i cdk8.smi -d dbs/cdk8.dat -p my_project/ -c 4

# 3. Screen a new compound library with the generated models
screen_db -d dbs/new_library.dat -q my_project/models/ -o my_project/vs/ -c 4 -v
```

The input file (`cdk8.smi`) is a tab-separated SMILES file with three columns:
`SMILES`, `compound_id`, `activity` (1 = active, 0 = inactive).

```
smiles                              mol_id          activity
O=C1NCCC12CCN(...)CC2               CHEMBL3798663   1
CC(C)(O)Cn1cc(...)cn1               CHEMBL3798944   1
c1ccc2ccccc2c1                      CHEMBL123456    0
```

Output pharmacophore models are plain `.xyz` files stored in
`my_project/models/`. Validation statistics are written to
`my_project/external_statistics.txt`.

## Project layout after a full run

```
my_project/
├── trainset/           Training set .smi files (centroids.smi, t0.smi, t1.smi, ...)
├── models/             Pharmacophore model .xyz files
├── raw_screen/         Hit lists used for internal validation
└── external_statistics.txt   Model performance metrics
```

## Cite

Ligand-Based Pharmacophore Modeling Using Novel 3D Pharmacophore Signatures
Kutlushina et al., *Molecules* **2018**, 23(12), 3094
https://doi.org/10.3390/molecules23123094

Probabilistic Approach for Virtual Screening Based on Multiple Pharmacophores
Madzhidov et al., *Molecules* **2020**, 25(2), 385
https://doi.org/10.3390/molecules25020385
