# mauve-py
Python wrappers for a fork of progressiveMauve multiple genome aligner

mauve-py is particularly suited to mapping modified synthetic genomes back
to a reference genome and includes features to aid in the "healing" of the
mapping that the progressive mauve binaries miss

See [mauve alignment](http://asap.ahabs.wisc.edu/software/mauve/overview.html) for original source. 
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
Copyright (c) 2014  Nick Conway, Ben Pruitt; Wyss Institute for 
Biologically Inspired Engineering


## Usage
Right now we support comparing only two genomes at a time, although mauve will
handle more.

Use fasta files on the input:

```python
import mauve
index_lookup_table = mauve.buildIndex("my_genome_A.fa", "reference_genome_B.fa")
```

the `index_lut` is a lookup table mapping a modified genome A index to a reference 
genome B index from the 5' end starting at 0.

## Dependencies

### Binary dependencies

Get these for your platform (Linux or OS X, Windows is unsupported)
    
- [automake](https://www.gnu.org/software/automake/)
- [libtool](https://www.gnu.org/software/libtool/)
- [pkg-config](http://www.freedesktop.org/wiki/Software/pkg-config/)
- [boost](http://www.boost.org/) > 1.55

Which can easily be installed with the package manager for your OS, 
apt, yum, etc on Linux, or [Homebrew](http://brew.sh/) for OS X

### Python dependencies

- [numpy](https://pypi.python.org/pypi/numpy/1.9.0)
- [bitarray](https://pypi.python.org/pypi/bitarray/0.8.1)
- [cython](https://pypi.python.org/pypi/Cython/0.21)
- [libnano](https://github.com/Wyss/libnano)

which can all be installed with `pip` except `libnano`

## LINUX INSTALLATION (Debian)
Double check you have the right `boost` components installed

`boost` requirements for debian linux

    libboost-filesystem-dev, libboost-program-options-dev, libboost-iostreams-dev

run:

    python setup.py install

## OS X INSTALLATION (Mavericks tested)
Double check you have the right boost components installed

make sure you have `Homebrew` installed on OS X Mavericks

    brew install boost --c++11

run: 

    python setup.py install

see:

-   [clang boost ref 1](http://hnrkptrsn.github.io/2013/02/26/c11-and-boost-setup-guide/)
-   [clang boost ref 2](http://stackoverflow.com/questions/17884344/why-does-boost-compilation-fails-with-clang)

# Building Mauve stand-alone

If you'd like to use mauve without the python wrapper, look in the `mauve/src`
path for the fork of the source.  It is updated to build on recent gcc and clang
with c++11, mainly out of the need to get rid of build errors and warnings.

## Deviations from stock mauve
    
- updated to latest boost 1.55ish so we now run BOOST_FILESYSTEM_VERSION 3 instead of 2
- OS X is clang compatible now. follow instructions below.
- no Windows support although it might be easy to add back in  
- lots of compiler warning fixes thanks to building on both clang and gcc