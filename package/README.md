Install Models with:
```
nucleofind-install -m phos -o site_packages
```

Example Usage:

```
nucleofind -i PDB.mtz -o PDB.map -intensity FWT -phase PHWT
```


```
usage: nucleofind [-h] [-m M] -i I -o O [-r [R]] [-intensity [INTENSITY]] [-phase [PHASE]]

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
usage: nucleofind-install [-h] -m {phos,sugar,base} [-o {site_packages,ccp4}] [--all] [--reinstall]

nucleofind Install

options:
  -h, --help            show this help message and exit
  -m {phos,sugar,base}, --model {phos,sugar,base}
  -o {site_packages,ccp4}, --output {site_packages,ccp4}
  --all
  --reinstall
```