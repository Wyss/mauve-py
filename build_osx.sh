#!/bin/sh

# brew install automake
# brew install boost --c++11
# brew install pkg-config
# brew install libtool

export CFLAGS="-O3 -arch x86_64 -mmacosx-version-min=10.9"
export CXXFLAGS="-O3 -arch x86_64 -mmacosx-version-min=10.9 -DBOOST_SYSTEM_NO_DEPRECATED "
export CXX="clang++ -std=c++11 -stdlib=libc++"
# export BJAMFLAGS=" --toolset=darwin --user-config=${BPATH}/user-config.jam  link=static threading=single "
BPATH=`pwd`
export BJAMFLAGS=" --toolset=darwin --user-config=${BPATH}/user-config.jam link=static"
echo $BJAMFLAGS
export CONFIGFLAGS=" --disable-shared --disable-dependency-tracking "
export PATH=$PATH:"/usr/local/lib:/usr/local/include"
# export PATCHBOOST="patch"
# export PATCHARGS=" -p1 -i ../../details/boost_1_38_0-mac_os_x_gcc42.patch"
export ANTDISTCMD="ant macdist"

./build_mauve.sh
