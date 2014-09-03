See [mauveAligner](http://asap.ahabs.wisc.edu/software/mauve/overview.html) for original source. 
The Mauve License is GPLv2.  

Please cite for progressive mauve:

    Aaron E. Darling, Bob Mau, and Nicole T. Perna. 2010.  progressiveMauve: Multiple Genome Alignment with Gene Gain, Loss, and Rearrangement.  PLoS One.  5(6):e11147. 

All dependencies are grouped in this repo which include [muscle](http://www.drive5.com/muscle/).
Muscle is licensed to the [Public Domain](http://www.drive5.com/muscle/manual/license.html)

For Muscle cite:

    Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy 
    and high throughput Nucleic Acids Res. 32(5):1792-1797

The included Python wrappers and file parsers are GPLv2 Licensed
Copyright (c) 2014 Ben Pruitt, Nick Conway; Wyss Institute for 
Biologically Inspired Engineering

# Deviations from stock mauve
    
updated to latest boost 1.55ish so we now run BOOST_FILESYSTEM_VERSION 3 instead of 2

os x is clang compatible now. follow instructions below  

# dependencies

Get these for your platform (Linux or OS X)
    
    automake
    libtool
    pkg-config
    boost > 1.55ish

# INSTALLATION preliminary

set the MAUVE_DIR environment variable to path to this repository on the 
installation machine

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
associated with missing alignements in module buildindex which depends on 

- numpy
- bitarray
- cython 

run:

    python setup.py build_ext --inplace

to build inplace or install them with:

    python setup.py install