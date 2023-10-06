from setuptools import setup
import os
import argparse

def read_file(filename):
    with open(os.path.join(os.path.dirname(__file__), filename)) as file:
        return file.read()

setup(
    name="xtal-cartographer",
    author="Jordan Dialpuri",
    author_email="jsd523@york.ac.uk",
    version="0.1.6",
    entry_points={
        "console_scripts": [
            'cartographer=cartographer.predict:run',
            'cartographer-install=cartographer.install:run'
        ]
    },
     install_requires=[ 
        'tensorflow',
        'tensorflow_addons',
        'gemmi',
        'tqdm' ],
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',

)
