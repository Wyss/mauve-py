#!/bin/sh

# see http://asap.ahabs.wisc.edu/mauve/mauve-developer-guide/compiling-mauvealigner-from-source.html
# need boost > 1.34 and pkg-config
# curl http://voxel.dl.sourceforge.net/project/boost/boost/1.38.0/boost_1_38_0.tar.bz2 > boost_1_38_0.tar.bz2 2> /dev/null

# If you ever happen to want to link against installed libraries
# in a given directory, LIBDIR, you must either use libtool, and
# specify the full pathname of the library, or use the `-LLIBDIR'
# flag during linking and do at least one of the following:
#    - add LIBDIR to the `LD_LIBRARY_PATH' environment variable
#      during execution
#    - add LIBDIR to the `LD_RUN_PATH' environment variable
#      during linking
#    - use the `-Wl,-rpath -Wl,LIBDIR' linker flag
#    - have your system administrator add LIBDIR to `/etc/ld.so.conf'

VERSION="2-1-0"
BPATH=`pwd`
export PATH=$PATH:$BPATH/bin
export PKG_CONFIG_PATH=$BPATH/lib/pkgconfig


# svn co https://svn.code.sf.net/p/mauve/code/libGenome/trunk libGenome
# svn co https://svn.code.sf.net/p/mauve/code/muscle/trunk muscle
# svn co https://svn.code.sf.net/p/mauve/code/libMems/trunk libMems
# svn co https://svn.code.sf.net/p/mauve/code/mauveAligner/trunk mauveAligner

cd libGenome
./autogen.sh && ./configure --prefix=$BPATH && make clean && make -j2 install
cd ..


cd muscle
./autogen.sh && ./configure --prefix=$BPATH && make clean && make -j2 ; make && make install
cd ..

cd libMems
# ./autogen.sh && ./configure --prefix=$BPATH --with-boost=$HOME && make clean && make -j2 install
./autogen.sh && ./configure --prefix=$BPATH && make clean && make -j2 install
cd ..

cd mauveAligner
./autogen.sh && ./configure --prefix=$BPATH && cd src && make clean && make mauveStatic && make progressiveMauveStatic
cd ..


