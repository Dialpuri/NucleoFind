# Cartographer

[![PyPI version](https://badge.fury.io/py/xtal-cartographer.svg)](https://badge.fury.io/py/xtal-cartographer)

Nucleic acid electron density interpretation remains a difficult problem for computer programs to deal with. Programs tend to rely on exhaustive searches to recognise characteristic features. Cartographer is a deep-learning-based approach to interpreting and segmenting electron density. Using a crystallographic map, the positions of the phosphate group, sugar ring and nitrogenous base group are able to be predicted with high accuracy. 

## Requirements
- python3

## User Installation 
```
pip install xtal-cartographer
cartographer-install -o site_packages --all
```

## Usage

```
cartographer -i PDB.mtz -o PDB.map -intensity FWT -phase PHWT
```


```
usage: cartographer [-h] [-m M] -i I -o O [-r [R]] [-intensity [INTENSITY]] [-phase [PHASE]]

options:
  -h, --help              show this help message and exit
  -m M, -model_path M     Path to model
  -i I, -input I          Input mtz path
  -o O, -output O         Output map path 
  -r [R], -resolution [R] Resolution cutoff to apply to mtz
  -intensity [INTENSITY]  Name of intensity column in MTZ
  -phase [PHASE]          Name of phase column in MTZ
```

```
usage: cartographer-install [-h] -m {phos,sugar,base} [-o {site_packages,ccp4}] [--all] [--reinstall]

Cartographer Install

options:
  -h, --help            show this help message and exit
  -m {phos,sugar,base}, --model {phos,sugar,base}
  -o {site_packages,ccp4}, --output {site_packages,ccp4}
  --all
  --reinstall
```

## Developer Installation 
Clone the project

```
git clone https://github.com/Dialpuri/Cartographer.git
```

Change directories into Cartographer

```
cd Cartographer
```

Create a Python virtual environment and entire the environment

```
python3 -m virtualenv pyenv
source pyenv/bin/activate
```
Install using pip

```
cd cartographer
pip install .
```

