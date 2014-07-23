#!/bin/sh
mkdir -p config
autoreconf --force --install -I config  
echo "Now run ./configure --prefix=$HOME ; make install"

