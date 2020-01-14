# PSearch - 3D ligand-based pharmacophore modeling

PSearch is a tool to generate 3D ligand-based pharmacophore models and perform virtual screening with them.

## Installation

```bash
pip install psearch
```

## Dependency

`pmapper >= 0.3.1`

## Example

### Creation of ligand-based pharmacophore models
It is recommended to create an empty dir which would be your `$PROJECT_DIR` and copy an input file to that location.  
There are two steps of pharmacophore model generation.  

1. Dataset preparation. 

```python
prepare_datatset -i $PROJECT_DIR/input.smi -c 4
```
`-i` - path to the input file;  
`-c` - number of CPUs to use.
There are some other arguments which one can use. Invoke script with `-h` key to get full information.  

The script takes as input a tab-separated SMILES file containing `SMILES`, `compound id`, `activity` columns without a header. 
The third column should contain a word `active` or `inactive`.
The script splits input compounds on active and inactive subsets, generates stereoisomers and conformers, creates databases of active and inactive compounds with labeled pharmacophore features.  

2. Model building.  

```python
psearch -p $PROJECT_DIR -t 0.4 -ts 1 2 -c 4
```
`-p` - path to the project dir;  
`-t` - threshold for compound clustering to create training sets;
`-ts` - modes of formed training sets, 1 - form a training set by Strategy 1 (a single training set from centroids of individual clusters), 2 - form a training set by Strategy 2 (separate training set per each cluster), 1 2 - form a training sets by Strategy 1 and Strategy 2;
`-c`- number of CPUs to use

### Virtual screening with pharmacophore models 

1. Database creation. 

The script takes as input a tab-separated SMILES file containing `SMILES` and `compound id` columns.

```python
prepare_db -i compounds.smi -o compounds.db -c 4 -v
```
`-i` - path to the input file;  
`-c` - number of CPUs to use;
`-v` - print progress 
There are other arguments available to tweak database generation. To get the full list of arguments invoke `-h` key.
 
2. Virtual screening.
  
```python
screen_db -d compounds.db -q $PROJECT_DIR/models/ -o screen_results/ -c 4
```
`-d` - input generated SQLite database  
`-q` - pharmacophore model or models or a directory with models   
If a directory would be specified all pma- and xyz-files will be recognized as pharmacophores and will be used for screening  
`-o` - path to an output directory if multiple models were supplied for screening or a path to a text file    
`-c`- number of CPUs to use

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
