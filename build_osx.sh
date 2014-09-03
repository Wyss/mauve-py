#!/bin/sh

# brew install automake
# brew install boost --c++11
# brew install pkg-config
# brew install libtool

export CFLAGS="-O3 -arch x86_64 -mmacosx-version-min=10.9 "
export CXXFLAGS="-O3 -arch x86_64 -mmacosx-version-min=10.9 -fcxx-exceptions -DBOOST_SYSTEM_NO_DEPRECATED -DHAVE_MKSTEMP "
export CXX="clang++ -std=c++11 -stdlib=libc++"

BPATH=`pwd`
export BJAMFLAGS=" --toolset=darwin --user-config=${BPATH}/user-config.jam link=static"
echo $BJAMFLAGS
export CONFIGFLAGS=" --disable-shared --disable-dependency-tracking "
export PATH=$PATH:"/usr/local/lib:/usr/local/include"

export ANTDISTCMD="ant macdist"

./build_mauve.sh
