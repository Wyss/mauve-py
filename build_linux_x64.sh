#!/bin/sh

# boost requirements for debian linux
# libboost-filesystem-dev, libboost-program-options-dev, libboost-iostreams-dev
# works when building with
# ./autogen.sh
# ./configure --prefix=$MAUVE_HOME
# make -j2
# make install

# build order:
# 1. libGenome
# 2. muscle
# 3. libMems
# 4. mauveAligner
which bjam > /dev/null
if [ $? -ne 0 ]; then 
  "Error: Boost's bjam must be installed to run this script\n"; exit;
  exit 1;
fi
which pkg-config > /dev/null
if [ $? -ne 0 ]; then 
  "Error: pkg-config must be installed to run this script\n"; exit;
  exit 1;
fi


mkdir -p build && \
cd build && \
export NO_COMPRESSION="true" && \
export ANTDISTCMD="ant dist" && \
../release_build.sh && \
cd .. && \
echo "done building binaries!\n"

