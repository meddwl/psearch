# PSearch - 3D ligand-based pharmacophore modeling

PSearch is a tool to automatically generate 3D ligand-based pharmacophore models.

## Installation

```bash
pip install psearch
```

## Dependency

`pmapper >= 0.2`  

## Example

### Creation of ligand-based pharmacophore models
It is recommended to create an empty dir which would be your `$PROJECT_DIR` and copy an input file to that location.  
There are two steps of pharmacophore model generation.  

1. Data set preparation. It takes as input a comma-separated SMILES file containing `SMILES`, `compound id`, `activity value`. It splits the input on active and inactive subsets, generates stereoisomers and conformers, creates databases of active and inactive compounds with labeled pharmacophore features.
```python
python3 prepare_datatset.py -i $PROJECT_DIR/input.smi -l 6 -u 8 -c 4
```
`-i` - path to the input file;  
`-u` - threshold to define active compounds (compounds with `activity value >= threshold` are considered active);  
`-l` - threshold to define inactive compounds (compounds with `activity value <= threshold` are considered inactive);  
`-c` - number of CPUs to use.  
There are other arguments available to tweak data set preparation. To get the full list of arguments run `python3 prepare_datatset.py -h`  

2. Model building.  

```python
python3 psearch.py -p $PROJECT_DIR -t 0.4 -c 4
```
`-p` - path to the project dir;  
`-t` - threshold for compound clustering to create training sets;  
`-c`- number of CPUs to use

### Virtual screening with pharmacophore models 

TODO

## Documentation

All scripts have `-h' argument to retrieve descriptions of all available options and arguments.

## Authors
Alina Kutlushina, Pavel Polishchuk

## Citation
Ligand-Based Pharmacophore Modeling Using Novel 3D Pharmacophore Signatures  
Alina Kutlushina, Aigul Khakimova, Timur Madzhidov, Pavel Polishchuk  
*Molecules* **2018**, 23(12), 3094  
https://doi.org/10.3390/molecules23123094

## License
BSD-3 clause
