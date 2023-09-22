from setuptools import setup

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
     install_requires=[
        'tensorflow',
        'tensorflow_addons',
        'gemmi',
        'tqdm'
    ]
)
