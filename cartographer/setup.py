from setuptools import setup
import urllib.request
import site 
import os

def download_dependencies(): 
    site_packages_dir = site.getsitepackages()
    if site_packages_dir:
        model_dir = os.path.join(site_packages_dir[0], "cartographer_models")
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)

        phos_model_path = os.path.join(model_dir, "phos.hdf5")
        if not os.path.exists(phos_model_path):
            urllib.request.urlretrieve("https://github.com/Dialpuri/Cartographer/raw/master/models/phos.hdf5", phos_model_path)  

        base_model_path = os.path.join(model_dir, "base.hdf5")
        if not os.path.exists(base_model_path):
            urllib.request.urlretrieve("https://github.com/Dialpuri/Cartographer/raw/master/models/base.hdf5", base_model_path)

    return [ 
        'tensorflow',
        'tensorflow_addons',
        'gemmi',
        'tqdm', ]

setup(
    name="Cartographer",
    author="Jordan Dialpuri",
    author_email="jsd523@york.ac.uk",
    version="0.1",
    entry_points={
        "console_scripts": [
            'cartographer=cartographer.predict:run'
        ]
    },
     install_requires=download_dependencies(),

)
