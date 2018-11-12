# PSearch - 3D ligand-based pharmacophore modeling

PSearch is a Python application to automatically generate 3D pharmacophore models based on a supplied data set of compounds with measured activity values.

## Installation

```bash
git clone https://github.com/meddwl/psearch
git submodule init
git submodule update
```

## Dependency

`rdkit >= 2017.09`  
`networkx >= 1.11`  

## Example

It is recommended to create an empty dir which would be your `$PROJECT_DIR` and copy an input file to that location.  
There are two steps of pharmacophore model generation.  

'I.' 

1. Data set preparation. It takes as input a comma-separated SMILES file containing `SMILES`, `compound id`, `activity value`. It splits the input on active and inactive subsets, generates stereoisomers and conformers, creates databases of active and inactive compounds with labeled pharmacophore features.
```python
python3 prepare_datatset.py -i $PROJECT_DIR/input.smi -l 6 -u 8 -c 4
```
`-i` - path to the input file;  
`-u` - treshold to define active compounds (compounds with `activity value >= threshold` are considered active);  
`-l` - treshold to define inactive compounds (compounds with `activity value <= threshold` are considered inactive);  
`-c` - number of CPUs to use.  
There are other arguments available to tweak data set preparation. To get the full list of agruments run `python3 prepare_datatset.py -h`  

2. Model building.  

```python
python3 psearch.py -p $PROJECT_DIR -t 0.4 -c 4
```
`-p` - path to the project dir;  
`-t` - threshold for compound clustering to create training sets;  
`-c`- number of CPUs to use

'II.' 

## Authors
Alina Kutlushina, Pavel Polishchuk

## Citation
...

## License
BSD-3
