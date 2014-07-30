#!/bin/sh

# see http://asap.ahabs.wisc.edu/mauve/mauve-developer-guide/compiling-mauvealigner-from-source.html
# need boost > 1.34 and pkg-config
# curl http://voxel.dl.sourceforge.net/project/boost/boost/1.38.0/boost_1_38_0.tar.bz2 > boost_1_38_0.tar.bz2 2> /dev/null

VERSION="2-1-0"
BPATH=`pwd`
export PATH=$PATH:$BPATH/bin
export PKG_CONFIG_PATH=$BPATH/lib/pkgconfig


# svn co https://svn.code.sf.net/p/mauve/code/libGenome/trunk libGenome
# svn co https://svn.code.sf.net/p/mauve/code/muscle/trunk muscle
# svn co https://svn.code.sf.net/p/mauve/code/libMems/trunk libMems
# svn co https://svn.code.sf.net/p/mauve/code/mauveAligner/trunk mauveAligner

# build order:
# 1. libGenome
# 2. muscle
# 3. libMems
# 4. mauveAligner

echo "Begin libGenome"
cd libGenome
./autogen.sh && ./configure --prefix=$BPATH && make clean && make -j2 install
cd ..
echo "End libGenome"

echo "Begin muscle"
cd muscle
./autogen.sh && ./configure --prefix=$BPATH && make clean && make -j4 ; make && make install
cd ..
echo "End muscle"

echo "Begin libMems"
cd libMems
# ./autogen.sh && ./configure --prefix=$BPATH --with-boost=$HOME && make clean && make -j2 install
# ./autogen.sh && ./configure --prefix=$BPATH && make clean && make -j2 install
./autogen.sh && ./configure --prefix=$BPATH && make clean && make -j4 install
cd ..
echo "End libMems"

echo "Begin mauveAligner"
cd mauveAligner
# ./autogen.sh && ./configure --prefix=$BPATH && cd src && make clean && make progressiveMauveStatic install
./autogen.sh && ./configure --prefix=$BPATH && cd src && make clean && make -j4 progressiveMauveStatic install
echo "End mauveAligner"

