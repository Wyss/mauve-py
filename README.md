see [mauveAligner](http://asap.ahabs.wisc.edu/mauve/mauve-developer-guide/compiling-mauvealigner-from-source.html) for original source. License is GPLV2.  
All dependencies are grouped in this repo

# Deviations from stock mauve
    
updated to latest boost 1.55ish so we now run BOOST_FILESYSTEM_VERSION 3 instead of 2

os x is clang compatible now. follow instructions below  

# dependencies

Get these for your platform (Linux or OS X)
    
    automake
    libtool
    pkg-config
    boost > 1.55ish

# LINUX installation (Debian)

boost requirements for debian linux

    libboost-filesystem-dev, libboost-program-options-dev, libboost-iostreams-dev

run

    ./build_mauve.sh

# OS X INSTALLATION
make sure you have Homebrew installed on OS X Mavericks

    brew install boost --c++11

run  

    ./build_osx.sh

see:

-   [clang boost ref 1](http://hnrkptrsn.github.io/2013/02/26/c11-and-boost-setup-guide/)
-   [clang boost ref 1](http://stackoverflow.com/questions/17884344/why-does-boost-compilation-fails-with-clang)