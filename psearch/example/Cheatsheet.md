# PSearch - 3D ligand-based pharmacophore modeling

PSearch is a tool to generate 3D ligand-based pharmacophore models and perform virtual screening with them.

## Installation

```bash
pip install -U git+https://github.com/meddwl/psearch.git@dev2
```

## Dependency

`pmapper >= 0.4.0`

## Example

1. Generation of a database with precomputed conformers and pharmacophores.

```python
gen_db -i Z2044051982.smi -o data/Z2044051982.dat -n 50 -c 4
```
`-i` - path to the input SMILES file
`-o` - path to database (should have extension .dat) 
`-n` - number of generated conformers (the minimum recommended number of generated conformers is 50)
`-c` - number of CPUs to use  

The script takes as input a tab-separated SMILES file containing `SMILES`, `compound id` columns. 
The script generates stereoisomers and conformers, creates the database of compounds with labeled pharmacophore

2. Virtual screening.
  
```python
screen_db -d data/Z2044051982.dat -q psearch/chembl_models/CHEMBL204 psearch/chembl_models/CHEMBL205 psearch/chembl_models/CHEMBL3975 -o screen_results/ -c 4
```
`-d` - input generated database  
`-q` - pharmacophore model or models or a directory with models. If a directory would be specified all pma- and xyz-files will be recognized as pharmacophores and will be used for screening.  
`-o` - path to an output directory if multiple models were supplied for screening or a path to a text/sdf file    
`-c`- number of CPUs to use

If multiple models are used for screening and sdf output is desired a user should add `-z` argument which will force output format to be sdf.

3. Calculate the predicted probability of the molecules

```python
prediction -vs screen_results/ -p pharmacophores_stat.csv -o resilts.txt
```
`-vs` - path to the virtual screening result  
`-q` - .csv file with the precision of pharmacophore models
`-o` - output text file where will be saved the prediction
