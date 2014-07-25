#!/bin/sh

# automake
# boost
# pkg-config
# libtool

export CFLAGS="-O3"
export CXXFLAGS="-O3 -DBOOST_SYSTEM_NO_DEPRECATED "
export BJAMFLAGS=" link=static "
echo $BJAMFLAGS
export CONFIGFLAGS=" --disable-shared --disable-dependency-tracking "
export PATH=$PATH:"/usr/lib:/usr/include"

./build_mauve.sh
