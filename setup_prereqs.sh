#!/bin/sh

source ./env.sh

cd eigen
rm -rf build
mkdir build
cd build
cmake -DCMAKE_CXX_FLAGS="-std=c++17 " -DCMAKE_INSTALL_PREFIX=$(pwd)/../.. ..
make install 

cd ../..
cd ifopt 
mkdir build
cd build
cmake -DCMAKE_CXX_FLAGS="-std=c++17 " -DCMAKE_INSTALL_PREFIX=$(pwd)/../.. ..
make install
cd ../..

