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
            print("Downloading the Cartographer phosphate model into ", phos_model_path, " cancel now if you do not wish for this to happen")
            urllib.request.urlretrieve("https://github.com/Dialpuri/Cartographer/raw/master/models/phos.hdf5", phos_model_path)  

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
