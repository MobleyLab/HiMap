"""
Setup script for HiMap

You can install himap with

python setup.py install
"""

import sys,os
from os.path import relpath, join
import versioneer

from setuptools import setup, find_packages

if sys.argv[-1] == 'setup.py':
    print("To install, run 'python setup.py install'")
    print()


descr = """
High Information Mapper (HiMap), successor of the Lead 
Optimization Mapper (LOMAP) is an automated algorithm
to plan efficient relative free energy calculations between
potential ligands within a substantial of compounds'
"""

setup(
    name                 = 'himap',
    version              = versioneer.get_version(),
    cmdclass             = versioneer.get_cmdclass(),
    description          = 'High Information Mapper',
    long_description     = descr,
    classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics'
    ],
    keywords=[ 'alchemical free energy setup', 'perturbation network' ],
    url                  = 'https://github.com/MobleyLab/HiMap',
    author               = 'Mary Pitman, Gaetano Calabro, David Mobley',
    maintainer           = 'Mary Pitman and David Mobley',
    author_email         = 'mpitman -at- uci.edu',
    license              = 'GNU 3.0',
    platforms            = ['Linux-64', 'Mac OSX-64', 'Unix-64'],
    packages             = find_packages()+['test'],
    include_package_data = True,

    entry_points         = {'console_scripts':['lomap=lomap.dbmol:startup']},
    zip_safe             = False
)

