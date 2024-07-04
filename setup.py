# setup.py

from setuptools import setup, find_packages

setup(
    name='band_lithium',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'ase',
        'pymatgen'
    ],
    entry_points={
        'console_scripts': [
            'band_lithium=band_lithium.core:main',
        ],
    },
)

