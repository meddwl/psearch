# Installation

## System requirements

- **Python** >= 3.9
- **Operating system**: platform-independent (Linux, macOS, Windows)
- RDKit is the largest dependency; it is most reliably installed through conda

## Recommended: conda environment

RDKit does not have an official PyPI wheel, so a conda environment is the
most reliable way to satisfy all dependencies.

```bash
# Create a new environment (Python 3.10 is a safe choice)
conda create -n psearch python=3.10

# Activate it
conda activate psearch

# Install RDKit from conda-forge
conda install -c conda-forge rdkit

# Install networkx (graph library used by pmapper internally)
conda install -c conda-forge networkx
```

## Install PSearch

**From PyPI** (stable release):

```bash
pip install psearch
```

**From GitHub** (latest development version):

```bash
pip install -U git+https://github.com/meddwl/psearch.git
```

Both commands pull in `pmapper >= 0.4.1` automatically as the only declared
Python dependency. RDKit and networkx must be installed separately via conda
(see above).

## Full dependency list

| Package | Minimum version | How to install |
|---------|----------------|---------------|
| Python | 3.9 | conda |
| rdkit | 2017.09 | `conda install -c conda-forge rdkit` |
| networkx | 2 | `conda install -c conda-forge networkx` |
| pmapper | 0.4.1 | installed automatically by pip |
| pandas | any recent | installed automatically |
| numpy | any recent | installed automatically |

## Verify the installation

Check that the CLI entry points are available:

```bash
gen_db -h
psearch -h
screen_db -h
prediction -h
external_stat -h
```

Or verify from Python:

```python
from psearch.database import DB
from psearch.gen_db import create_db
print("psearch imported successfully")
```

## Common installation issues

**`rdkit` not found after `pip install psearch`**
RDKit cannot be installed via pip alone. Install it through conda-forge first,
then run `pip install psearch` inside the same activated conda environment.

**`ModuleNotFoundError: No module named 'pmapper'`**
pmapper should be installed automatically by pip. If it is missing, run
`pip install pmapper`.

**Database files are not portable between platforms**
PSearch databases use Python's `shelve` module, which creates `.dat` and `.dir`
sidecar files. These files may not be readable across different operating
systems or Python versions. Build the database on the platform you intend to
use for screening.

**Database already exists error**
`gen_db` refuses to overwrite an existing database. Either choose a new output
path with `-d` or remove the existing `.dat` and `.dir` files manually before
re-running.
