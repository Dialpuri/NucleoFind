# NucleoFind <img src="https://github.com/Dialpuri/NucleoFind/assets/44945647/a7c6c30c-a9fb-4f30-bdc0-3705ae9df36f" alt="logo" width="50"/> 

[![PyPI version](https://badge.fury.io/py/nucleofind.svg)](https://badge.fury.io/py/nucleofind)
![GitHub Issues](https://img.shields.io/github/issues-raw/Dialpuri/NucleoFind)
![PyPI - Downloads](https://img.shields.io/pypi/dm/NucleoFind)
![GitHub repo size](https://img.shields.io/github/repo-size/Dialpuri/NucleoFind)
[![Build documentation](https://github.com/Dialpuri/NucleoFind/actions/workflows/deploy.yml/badge.svg)](https://github.com/Dialpuri/NucleoFind/actions/workflows/deploy.yml)

Nucleic acid electron density interpretation remains a difficult problem for computer programs to deal with. Programs tend to rely on exhaustive searches to recognise characteristic features. NucleoFind is a deep-learning-based approach to interpreting and segmenting electron density. Using a crystallographic map, the positions of the phosphate group, sugar ring and nitrogenous base group are able to be predicted with high accuracy. 

## Documentation
[Link to Documentation](https://dialpuri.github.io/NucleoFind/about-nucleofind.html)

## Development

To build NucleoFind upon every Python import, which will be useful in development, run: 

    pip install --no-build-isolation --config-settings=editable.rebuild=true -Cbuild-dir=build -ve .