#!/bin/sh
export CFLAGS="-O3 -arch ppc -arch ppc64 -arch i386 -arch x86_64 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -mmacosx-version-min=10.4"
export CXXFLAGS="-O3 -arch ppc -arch ppc64 -arch i386 -arch x86_64 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -mmacosx-version-min=10.4"
export PKG_CONFIG_PATH="`pwd`/osx_build/lib/pkgconfig"
export BJAMFLAGS=" --toolset=darwin --user-config=../user-config.jam link=static threading=single " && \
export CONFIGFLAGS=" --disable-shared --disable-dependency-tracking "
export PATCHBOOST="patch"
export PATCHARGS=" -p1 -i ../../details/boost_1_38_0-mac_os_x_gcc42.patch"
export ANTDISTCMD="ant macdist"
mkdir -p build && \
cd build && \
echo "using darwin : 4.2 : g++ : <cxxflags>\"-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch ppc -arch ppc64 -arch i386 -arch x86_64 -mmacosx-version-min=10.4\" <cflags>\"-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch ppc -arch ppc64 -arch i386 -arch x86_64 -mmacosx-version-min=10.4\" ;" >user-config.jam && \
../release_build.sh && \
cd .. 

