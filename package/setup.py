from setuptools import setup
import os
import argparse

def read_file(filename):
    with open(os.path.join(os.path.dirname(__file__), filename)) as file:
        return file.read()

setup(
    name="nucleofind",
    description="NucleoFind: A Deep-Learning Network for Interpreting Nucleic Acid Electron Density",
    author="Jordan Dialpuri",
    author_email="jsd523@york.ac.uk",
    version='0.4.7',
    entry_points={
        "console_scripts": [
            'nucleofind=nucleofind.predict:run',
            'nucleofind-install=nucleofind.install:run',
            'nucleofind-clean=nucleofind.clean:run'
        ]
    },
     install_requires=[ 
        'onnxruntime',
        'gemmi',
        'tqdm' ],
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    packages=['nucleofind']

)
