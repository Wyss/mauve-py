#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
===========================================================
 mauve-py: 
===========================================================

``mauve-py`` is a Python wrappers for a derivative of the
popular mauve genome aligner.

See README.md for more information for usage examples.

Installation
------------

``mauve-py`` has three dependencies and should compile on
most linux and OS X systems that are running Python 2.7, 3.3, or 3.4.

To build ``mauve-py`` within the package directory run::

  $ python setup.py build_ext --inplace

If you would like to install mauve-py in your local Python environment
you may do so using the setup.py script::

  $ python setup.py install

'''

import os
import stat
import platform
import subprocess

from distutils import log as dlog
from distutils.core import setup, Extension
from distutils.command import install_lib
from os.path import join as pjoin
from os.path import relpath as rpath

dlog.set_verbosity(dlog.INFO)

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

try:
    import bitarray
except:
    print("Please install bitarray")
    raise


# list of all executable built by the mauve source
mauve_executables = ["addUnalignedIntervals", "makeBadgerMatrix",
"alignmentProjector", "backbone_global_to_local",
"bbAnalyze", "bbFilter",           
"coordinateTranslate", "createBackboneMFA",
"extractBCITrees", "getAlignmentWindows",
"getOrthologList", "repeatoire",
"mauveAligner", "mauveStatic",
"mauveToXMFA", "mfa2xmfa",
"muscle", "progressiveMauve",
"progressiveMauveStatic",
"projectAndStrip", "randomGeneSample",
"scoreAlignment", "stripGapColumns",
"stripSubsetLCBs", "toGrimmFormat",
"toMultiFastA", "toRawSequence",
"uniqueMerCount", "uniquifyTrees",
"xmfa2maf"
]


# to install executable files in package_data
class ExeInstallLib(install_lib.install_lib):
    def run(self):
        install_lib.install_lib.run(self)
        for fn in self.get_outputs():
            base_filename = os.path.basename(fn)
            if base_filename == "progressiveMauveStatic" or \
                    base_filename == "progressiveMauve":
                # make sure user executable
                mode = os.stat(fn).st_mode | stat.S_IXUSR | stat.S_IXGRP
                dlog.info("Changing mode of %s to %d", fn, mode)
                os.chmod(fn, mode)
# end class

with open('README.md') as fd:
    LONG_DESCRIPTION = fd.read()

PACKAGE_PATH =          os.path.abspath(os.path.dirname(__file__))
MODULE_PATH =           pjoin(PACKAGE_PATH, 'mauve')
SRC_PATH =              pjoin(MODULE_PATH, 'src')
BIN_PATH =              pjoin(SRC_PATH, 'bin')
INCLUDE_PATH =          pjoin(SRC_PATH, 'include')

# Build mauve for subprocess bindings

os_system = platform.system()
if os_system == 'Darwin':
    mauve_make = 'build_osx.sh'
elif os_system == 'Linux':      # update this for BSD, other POSIX etc
    mauve_make = 'build_linux.sh'
else:
    raise OSError("Platform %s not supported" % (os_system))

print("Building mauve dependencies...")
mauve_make = pjoin(SRC_PATH, mauve_make)
mauvebuild = subprocess.Popen([mauve_make], shell=True, cwd=SRC_PATH)
mauvebuild.wait()
print("Finished building mauve dependencies")

# Find all mauve data files to include with the package
mauve_files = [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
                    os.walk(BIN_PATH) for f in files]
mauve_files += [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
                    os.walk(INCLUDE_PATH) for f in files]


DESCRIPTION = ("Python wrapper for progressive mauve genome aligner")

DISTNAME = 'mauve-py'
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

mauve_extensions = [Extension('mauve.indexutils', 
                                sources=['mauve/indexutils.pyx'])]
mauve_ext_list = cythonize(mauve_extensions)

setup(
    name=DISTNAME,
    maintainer=AUTHORS,
    packages=['mauve'],
    ext_modules=mauve_ext_list,
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    maintainer_email=EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    download_url=DOWNLOAD_URL,
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS,
    package_data={'mauve': mauve_files},
    cmdclass={'install_lib':ExeInstallLib}
)