# Virtual Screening

## Overview

Virtual screening (VS) in PSearch matches 3D conformers of database compounds
against pharmacophore models using a two-stage protocol:

1. **Fingerprint pre-filter** — the pharmacophore fingerprint of the query
   model must be a bitwise subset of the compound's stored fingerprint.
   Compounds that fail this check are skipped without 3D fitting (fast).
2. **3D pharmacophore fitting** — for compounds passing the pre-filter, every
   stored conformer is aligned to the model using pmapper's `fit_model`. A hit
   is recorded if the alignment succeeds within the fitting tolerance.

This two-stage design allows PSearch to screen large databases efficiently:
the fingerprint pre-filter typically eliminates the vast majority of compounds
before any geometry is tested.

---

## Preparing the screening database

A screening database must be built with `gen_db` before any VS run. The
`--bin_step` used here **must match** the bin step used when the pharmacophore
models were built. The default is 1 Å in both steps.

```bash
gen_db -i library.smi -d dbs/library.dat -c 4 -n 50 -v
```

The input SMILES file needs at least two columns (SMILES and mol_id); the
activity column is not required for screening-only workflows:

```
CC(=O)Oc1ccccc1C(=O)O    aspirin
c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34    pyrene
```

If the library contains duplicate SMILES or duplicate mol_ids, `gen_db` will
warn and write a corrected `*-updated.smi` alongside the original file. The
database consists of two sidecar files (`.dat` and `.dir`) that must be kept
together.

---

## Running virtual screening

### Screen with your own models

```bash
screen_db -d dbs/library.dat \
          -q my_project/models/ \
          -o my_project/vs/ \
          -c 4 -v
```

`-q` accepts:
- A single `.pma` or `.xyz` file → output is a single `.txt` hit-list file
- A directory → output is a directory of `.txt` files, one per model
- Multiple directories → each is mirrored as a subdirectory under the output path

### Screen with the bundled ChEMBL models (multi-target profiling)

Omit `-q` to use the pre-built ChEMBL pharmacophore library that ships with
PSearch:

```bash
screen_db -d dbs/library.dat -o profiling/vs/ -c 4 -v
```

---

## Key options

| Flag | Default | Purpose |
|------|---------|---------|
| `-z/--output_sdf` | off | Write matching 3D conformers to `.sdf` files alongside each hit list |
| `--conf` | off | Report every matching conformer separately, not just the first per compound (required for CCA scoring — see below) |
| `-f/--min_features` | None | Skip models with fewer than N distinct-coordinate features |
| `-c/--ncpu` | 1 | Parallel worker processes; scales well to all available cores |
| `-v/--verbose` | off | Print a progress line to stderr every 10 molecules |

---

## Hit list output format

Each hit list is a tab-separated text file with one line per hit:

```
mol_id          stereo_id   conf_id
CHEMBL3798663   0           12
CHEMBL3800311   0           7
```

- `stereo_id` — index of the enumerated stereoisomer (0-based)
- `conf_id` — index of the conformer within that stereoisomer (0-based)

Hit-list filenames follow the convention `<target_id>.<model_id>.txt`, which
is the format expected by the downstream `prediction` tool.

---

## SDF output

Add `-z` to write the matching 3D conformers (superimposed onto the
pharmacophore model) to an SDF file alongside each hit list. The SD property
`RMSD` records the root-mean-square deviation of the pharmacophore fit in Å:

```bash
screen_db -d dbs/library.dat -q my_project/models/ -o my_project/vs/ -z -c 4
```

This produces pairs of files such as:
```
my_project/vs/cdk8.t1_f7_p3.txt   — hit list
my_project/vs/cdk8.t1_f7_p3.sdf   — 3D structures of the hits
```

---

## Conformer-coverage approach (CCA)

By default `screen_db` stops at the first conformer of a compound that matches
a model (`match_first_conf=True`). The conformer-coverage approach instead
records **all** matching conformers, which can be used for alternative scoring
based on how many of a compound's conformers satisfy the pharmacophore. To
enable this mode, add `--conf`:

```bash
screen_db -d dbs/library.dat -q models/ -o vs/ --conf -c 4
```

With `--conf`, the same compound may appear multiple times in the hit list
(once per matching conformer). Downstream tools that compute a CCA score expect
this format.

---

## Filtering models before screening

If you have many models and want to focus on the most selective ones, use
`-f/--min_features` to skip models with fewer than N distinct-coordinate
features:

```bash
# Only screen models with 4 or more distinct features
screen_db -d dbs/library.dat -q models/ -o vs/ -f 4 -c 4
```

Alternatively, use `get_pharm_props` to compute feature counts for all models
and select a subset manually before running `screen_db`.

---

## Extracting matching 3D structures

Use `get_conf_and_pharm` to pull a specific conformer from the database as an
SDF file together with its pharmacophore representation in `.xyz` format. The
required IDs (mol_id, stereo_id, conf_id) come from the hit list.

```bash
get_conf_and_pharm -d dbs/library.dat \
                   -m CHEMBL3798663 \
                   -s 0 \
                   -c 12 \
                   -o extracted/
```

This writes:
- `extracted/CHEMBL3798663-s0-c12.sdf` — the 3D conformer
- `extracted/CHEMBL3798663-s0-c12` — pharmacophore features in XYZ format

Multiple compounds can be extracted in one call by supplying multiple values to
`-m`, `-s`, and `-c` (positionally matched):

```bash
get_conf_and_pharm -d dbs/library.dat \
                   -m CHEMBL3798663 CHEMBL3800311 \
                   -s 0 0 \
                   -c 12 7 \
                   -o extracted/
```

---

## Performance notes

- Use `-c` to parallelise across all available cores; screening is
  embarrassingly parallel at the compound level.
- The fingerprint pre-filter typically eliminates more than 90 % of compounds
  before any 3D fitting is attempted. Models with very few features (3) produce
  less selective fingerprints and therefore pass more compounds to the 3D
  fitting stage, making screening slower.
- Use `-f 4` or higher to skip very simple models when speed is critical.
- The database is opened read-only during screening; it is safe to run multiple
  `screen_db` processes against the same database file simultaneously.
- Pre-existing output files are **overwritten** if `screen_db` is run again
  with the same output path.
