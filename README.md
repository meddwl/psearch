# PSearch — 3D Ligand-Based Pharmacophore Modeling

PSearch is a tool for generating 3D ligand-based pharmacophore models and performing virtual screening. It enumerates stereoisomers, embeds 3D conformers, computes pharmacophore fingerprints, builds models from clustered active sets, and ranks screened molecules by predicted activity probability.

## Installation

```bash
# Latest release from PyPI
pip install psearch

# Development version from GitHub
pip install -U git+https://github.com/meddwl/psearch.git
```

## Requirements

- Python >= 3.10
- [RDKit](https://www.rdkit.org/)
- [pmapper](https://github.com/DrrDom/pmapper) >= 0.4.1

## Workflow Overview

```
Input SMILES
     │
     ▼
  gen_db          ← build conformer/pharmacophore database
     │
     ▼
  psearch         ← cluster actives, generate models, screen, validate
  (or separately:
    screen_db     ← screen a database with pharmacophore models
    external_stat ← compute external validation statistics
    prediction    ← rank molecules by predicted activity)
```

## Usage

### 1. Build a conformer/pharmacophore database — `gen_db`

Reads a 2D SMILES file, enumerates stereoisomers, embeds 3D conformers (ETKDGv3 + MMFF), and stores the results in a single-file SQLite database (`.db`).

```bash
gen_db -i molecules.smi -d dbs/molecules.db -c 4 -v
```

The input SMILES file must be **tab-separated** with three columns: `SMILES`, `compound_id`, `activity` (1 = active, 0 = inactive).

| Argument | Default | Description |
|---|---|---|
| `-i` / `--input` | required | Input 2D SDF or tab-separated SMILES file |
| `-d` / `--db` | required | Output database path (must have `.db` extension) |
| `-c` / `--ncpu` | 1 | Number of CPUs |
| `-n` / `--nconf` | 50 | Number of conformers per stereoisomer |
| `-s` / `--nstereo` | 5 | Maximum stereoisomers per compound |
| `-e` / `--energy_cutoff` | None | Discard conformers with MMFF energy > cutoff above lowest (kcal/mol) |
| `-r` / `--rms` | None | Discard conformers with pairwise RMS below cutoff (Å) |
| `-b` / `--bin_step` | 1 | Bin width (Å) for pharmacophore coordinate discretisation |
| `-p` / `--pharm_def` | None | Custom pmapper feature definition file |
| `--seed` | -1 | Random seed for conformer embedding (-1 = no seed) |
| `-v` / `--verbose` | False | Print progress to stdout |

If duplicate SMILES or compound IDs are detected in the input, a corrected `*-updated.smi` file is written alongside the original.

---

### 2. Build models and screen — `psearch`

Full pipeline: clusters active compounds, generates pharmacophore models, screens the database, and computes external validation statistics.

```bash
psearch -i molecules.smi -d dbs/molecules.db -p my_project/ -c 4
```

| Argument | Default | Description |
|---|---|---|
| `-i` / `--molecules` | required | Tab-separated SMILES file (SMILES, mol_id, activity) |
| `-d` / `--database` | required | Database built with `gen_db` |
| `-p` / `--project_dir` | auto | Project directory for all outputs |
| `-m` / `--mode_train_set` | [1, 2] | Training-set strategy: 1 = centroid set; 2 = one set per cluster |
| `-t` / `--threshold` | 0.4 | Butina clustering threshold |
| `-l` / `--lower` | 3 | Minimum number of pharmacophore features per model |
| `-u` / `--upper` | None | Maximum number of pharmacophore features per model |
| `-b` / `--bin_step` | 1 | Bin width (Å) |
| `-tol` / `--tolerance` | 0 | Tolerance for stereoconfiguration sign calculation |
| `--fcfp4` | False | Use FCFP4 fingerprints for clustering (default: pharmacophore FP) |
| `-c` / `--ncpu` | 1 | Number of CPUs |

**Project directory structure after a run:**

```
my_project/
├── trainset/               ← training set files per cluster
├── models/                 ← pharmacophore model files (.xyz)
│   └── <db>.t<n>_f<n>_p<n>.xyz
├── raw_screen/             ← per-model screening hit lists
└── external_statistics.txt ← validation metrics (precision, recall, …)
```

---

### 3. Screen a database — `screen_db`

Screens a database against one or more pharmacophore models. If no query is provided, the built-in ChEMBL pharmacophore models are used.

```bash
# Screen with custom models
screen_db -d dbs/molecules.db -q my_project/models/ -o my_project/vs/ -c 4 -v

# Screen with built-in ChEMBL models (multiprofiling)
screen_db -d dbs/molecules.db -o profiling/vs/ -c 4 -v
```

| Argument | Default | Description |
|---|---|---|
| `-d` / `--dbname` | required | Input database (`.db` file) |
| `-q` / `--query` | built-in | Model file(s) or directory; uses ChEMBL built-in models if omitted |
| `-o` / `--output` | required | Output file (`.txt`) or directory for hit lists |
| `-f` / `--min_features` | None | Skip models with fewer distinct-coordinate features than this |
| `-z` / `--output_sdf` | False | Write matching 3D conformers to SDF alongside hit lists |
| `--conf` | False | Report each matching conformer separately (required for CCA scoring) |
| `-c` / `--ncpu` | 1 | Number of CPUs |
| `-v` / `--verbose` | False | Print progress to stdout |

---

### 4. Predict activity — `prediction`

Computes the probability of activity for each molecule based on virtual screening results. Uses the precision of individual pharmacophore models for consensus scoring.

```bash
prediction -s my_project/vs/ -p my_project/external_statistics.txt -f mean -o results.txt
```

| Argument | Default | Description |
|---|---|---|
| `-s` / `--path_vs` | required | Directory with `screen_db` hit-list `.txt` files |
| `-p` / `--pharm_stat` | built-in | File with model precision statistics; uses built-in ChEMBL stats if omitted |
| `-f` / `--scoring_scheme` | `mean` | Consensus scoring: `max` = highest probability; `mean` = average |
| `-o` / `--output` | auto | Output TSV file for predictions |

---

### 5. Compute external validation statistics — `external_stat`

Calculates external validation metrics for a set of pharmacophore models against a labelled test set.

```bash
external_stat -i molecules.smi -t my_project/trainset/ -m my_project/models/ -s my_project/raw_screen/ -o stats.txt
```

| Argument | Default | Description |
|---|---|---|
| `-i` / `--molecules` | required | Tab-separated SMILES file (SMILES, mol_id, activity) |
| `-t` / `--trainset` | required | Directory with training set files |
| `-m` / `--models` | required | Directory with pharmacophore model files |
| `-s` / `--screen` | required | Directory with screening results |
| `-o` / `--output` | auto | Output TSV file for validation statistics |

---

## Example

The `example/` directory contains sample input files:

- `cdk8.smi` — 233 CDK8 ligands (CHEMBL5719) for model building
- `mols_for_profiling.smi` — molecules for multiprofiling against built-in ChEMBL models

**Ligand-based pharmacophore modeling:**

```bash
gen_db -i example/cdk8.smi -d dbs/cdk8.db -c 4 -v
psearch -i example/cdk8.smi -d dbs/cdk8.db -p my_project/ -c 4
```

**Multiprofiling against built-in ChEMBL models:**

```bash
gen_db -i example/mols_for_profiling.smi -d dbs/profiling.db -c 4 -v
screen_db -d dbs/profiling.db -o profiling/vs/ -c 4 -v
prediction -s profiling/vs/ -o profiling/results.txt
```

## Authors

Alina Denzler, Pavel Polishchuk

## Citation

Ligand-Based Pharmacophore Modeling Using Novel 3D Pharmacophore Signatures  
Alina Kutlushina, Aigul Khakimova, Timur Madzhidov, Pavel Polishchuk  
*Molecules* **2018**, 23(12), 3094  
https://doi.org/10.3390/molecules23123094

Probabilistic Approach for Virtual Screening Based on Multiple Pharmacophores  
Timur Madzhidov, Assima Rakhimbekova, Alina Kutlushina, Pavel Polishchuk  
*Molecules* **2020**, 25(2), 385  
https://doi.org/10.3390/molecules25020385

## License

BSD-3-Clause
