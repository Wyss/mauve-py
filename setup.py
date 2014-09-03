#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test with:
    python setup.py build_ext --inplace
"""
DESCRIPTION = ("**buildindex** is a Python package for aligning genomes using mauve aligner")
LONG_DESCRIPTION = """
**buildindex** is a Python package for aligning genomes using mauve aligner
License is GPLv2
"""

DISTNAME = 'buildindex'
LICENSE = 'GPLv2'
AUTHORS = "Nick Conway, Ben Pruitt"
EMAIL = "nick.conway@wyss.harvard.edu"
URL = ""
DOWNLOAD_URL = ''
CLASSIFIERS = [
    'Development Status :: 1 - Beta',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]

from distutils.core import setup, Extension
try:
    from Cython.Build import cythonize
except:
  print("Please install cython")
  raise

try:
    import numpy.distutils.misc_util
except:
    print("Please install numpy")
    raise

bi_extensions = [Extension('buildindex.indexutils', sources=['buildindex/indexutils.pyx'])]
bi_ext_list = cythonize(bi_extensions)

setup(
    name=DISTNAME,
    maintainer=AUTHORS,
    packages=['buildindex'],
    ext_modules=bi_ext_list,
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    maintainer_email=EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    download_url=DOWNLOAD_URL,
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS
)