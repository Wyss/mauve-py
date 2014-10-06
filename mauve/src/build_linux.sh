#!/bin/sh

# automake
# boost
# pkg-config
# libtool

export CFLAGS="-O3"
export CXXFLAGS="-O3 -std=c++11 -DBOOST_SYSTEM_NO_DEPRECATED  -DHAVE_MKSTEMP "
export BJAMFLAGS=" link=static "
echo $BJAMFLAGS
export CONFIGFLAGS=" --disable-shared --disable-dependency-tracking "
export PATH=$PATH:"/usr/lib:/usr/include"

./build_mauve.sh
