# Cartographer

## Requirements
The models used in this directory are hosted in GitHub large file storage, to download them as part of this directory you must have `git-lfs` installed. However, this is not required to use the package as they will be downloaded as part of the installation. 

## Installation
Clone the project

```
git clone https://github.com/Dialpuri/Cartographer.git
```

Change directories into Cartographer/cartographer

```
cd Cartographer/cartographer
```
Ensure you are in a python virtual environment and install using pip

```
pip install .
```
## Usage
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
