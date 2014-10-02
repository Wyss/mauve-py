See [mauveAligner](http://asap.ahabs.wisc.edu/software/mauve/overview.html) for original source. 
The Mauve License is GPLv2.  

Please cite progressive mauve:

    Aaron E. Darling, Bob Mau, and Nicole T. Perna. 2010.  
    progressiveMauve: Multiple Genome Alignment with Gene Gain, Loss, 
    and Rearrangement.  PLoS One.  5(6):e11147. 

All dependencies are grouped in this repo which include [MUSCLE](http://www.drive5.com/muscle/).
MUSCLE is licensed in the [Public Domain](http://www.drive5.com/muscle/manual/license.html)

For MUSCLE cite:

    Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy 
    and high throughput Nucleic Acids Res. 32(5):1792-1797

The included Python wrappers and file parsers are GPLv2 Licensed
Copyright (c) 2014 Ben Pruitt, Nick Conway; Wyss Institute for 
Biologically Inspired Engineering

# Deviations from stock mauve
    
- updated to latest boost 1.55ish so we now run BOOST_FILESYSTEM_VERSION 3 instead of 2
- os x is clang compatible now. follow instructions below.  
- lots of compiler warning fixes thanks to building on both clang and gcc

# dependencies

Get these for your platform (Linux or OS X, Windows is unsupported)
    
    automake
    libtool
    pkg-config
    boost > 1.55ish

# INSTALLATION preliminary

set the `MAUVE_DIR` environment variable to path to this repository on the 
installation machine.   The build scripts by default make the build install inplace.

# LINUX INSTALLATION (Debian)

boost requirements for debian linux

    libboost-filesystem-dev, libboost-program-options-dev, libboost-iostreams-dev

run

    ./build_linux.sh

# OS X INSTALLATION (Mavericks tested)
make sure you have Homebrew installed on OS X Mavericks

    brew install boost --c++11

run  

    ./build_osx.sh

see:

-   [clang boost ref 1](http://hnrkptrsn.github.io/2013/02/26/c11-and-boost-setup-guide/)
-   [clang boost ref 2](http://stackoverflow.com/questions/17884344/why-does-boost-compilation-fails-with-clang)

# Python installation

Included are Python wrappers for handling mauve aligner output and fixing errors
associated with gaps alignments where an actual alignment could be found.
Module buildindex does this which depends on:

- [numpy](https://pypi.python.org/pypi/numpy/1.9.0)
- [bitarray](https://pypi.python.org/pypi/bitarray/0.8.1)
- [cython](https://pypi.python.org/pypi/Cython/0.21)

which can all be installed with `pip`

run:

    python setup.py build_ext --inplace

to build inplace or install them with:

    python setup.py install


**Additionally** please remember to set the `MAUVE_DIR` environment 
variable to point to where you choose to.