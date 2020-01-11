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

1. Data set preparation. 

It takes as input a tab-separated SMILES file containing `SMILES`, `compound id`, `activity value` without header. It splits the input on active and inactive subsets, generates stereoisomers and conformers, creates databases of active and inactive compounds with labeled pharmacophore features.
```python
prepare_datatset -i $PROJECT_DIR/input.smi -c 4
```
`-i` - path to the input file;  
`-c` - number of CPUs to use.  
There are other arguments available to tweak data set preparation. To get the full list of arguments run `prepare_datatset -h`  

If you need to prepare a datatset and you don't need to split the input on active and inactive subsets you can use the command below. 
It takes as input a tab-separated SMILES file containing `SMILES`, `compound id`. It generates stereoisomers and conformers, creates databases of compounds with labeled pharmacophore features. 
```python
prepare_db -i $PROJECT_DIR/input.smi -c 4 -v
```
`-i` - path to the input file;  
`-c` - number of CPUs to use. 
`-v` - print progress 
There are other arguments available to tweak data set preparation. To get the full list of arguments run `prepare_datatset -h`  


2. Model building.  

```python
psearch -p $PROJECT_DIR -t 0.4 -ts 1 2 -c 4
```
`-p` - path to the project dir;  
`-t` - threshold for compound clustering to create training sets;
`-ts` - modes of formed training sets, 1 - form a training set by Strategy 1 (a single training set from centroids of individual clusters), 2 - form a training set by Strategy 2 (separate training set per each cluster), 1 2 - form a training sets by Strategy 1 and Strategy 2;
`-c`- number of CPUs to use
There are other arguments available to tweak data set preparation. To get the full list of arguments run `prepare_datatset -h`  

### Virtual screening with pharmacophore models 

1. Data set preparation. It takes as input a tab-separated SMILES file containing `SMILES`, `compound id`. It generates stereoisomers and conformers, creates databases of compounds with labeled pharmacophore features.

```python
prepare_db -i $PROJECT_DIR/input.smi -c 4 -v
```
`-i` - path to the input file;  
`-c` - number of CPUs to use;
`-v` - print progress 
There are other arguments available to tweak data set preparation. To get the full list of arguments run `prepare_datatset -h`  

2. Model building.  

```python
screen_db -d $PROJECT_DIR/databased.db -q $PROJECT_DIR/models/ -o $PROJECT_DIR/screen/ -c 4
```
`-d` - input SQLite database with generated conformers;  
`-q` - pharmacophore model or models or a directory path. If a directory is specified all pma- and xyz-files will be used for screening as pharmacophore models;
`-o` - path to an output directory or test (.txt);
`-c`- number of CPUs to use
There are other arguments available to tweak data set preparation. To get the full list of arguments run `screen_db -h`  

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
