# -*- coding: utf-8 -*-
from distutils.core import setup

setup(
    name='KiwiDist',
    version='0.3.6',
    author='F. Gatto and L. Väremo',
    author_email=['kiwi@sysbio.se'],
    packages=['kiwi'],
    license='LICENSE.txt',
    description='Combining gene-set analysis with network properties.',
    long_description=open('README.txt').read(),
    install_requires=[
        "matplotlib >= 1.3.1, <=1.4.3",
        "mygene >= 2.1.0",
        "networkx >= 1.8.1",
        "numpy >= 1.8.0",
        "pandas >= 0.13.1",
        "scipy >= 0.13.3, <=0.16.0",
        "six >=1.5"
        ],
)
