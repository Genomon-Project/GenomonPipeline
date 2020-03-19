#!/usr/bin/env python

from setuptools import setup, find_packages
from scripts.genomon_pipeline import __version__

#import sys
#sys.path.append('./tests')

setup(
    name = 'genomon_pipeline',
    version = __version__,
    description = 'Python tools for running genomon pipeline for cancer genome and transcriptome sequencing analysis',
    keywords = 'cloud bioinformatics',
    author = 'Kenichi Chiba, Ai Okada and Yuichi Shiraishi',
    author_email = 'genomon.devel@gmail.com',
    url = 'https://github.com/Genomon-Project/Genomon.git',
    license = 'License of GenomonPipeline',
    
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Genome Analysis :: RNA-seq',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    
    package_dir = {'': 'scripts'},
    packages = find_packages("scripts"),
    package_data = {'genomon_pipeline': ['*/data/*']},
    
    scripts = ['genomon_pipeline', 'genomon_runner'],
    include_package_data = True,
    zip_safe = False,
    install_requires = [
        'snakemake',
        'pyyml',
        'drmaa',
    ],
    #test_suite = 'unit_tests.suite'
    test_suite = 'tests'
)
