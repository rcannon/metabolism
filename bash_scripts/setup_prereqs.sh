#!/bin/bash

cd $PROJECT_DIR

mkdir PREREQS_INSTALL

echo
echo "Installing Eigen."
echo

cd eigen
rm -rf build
mkdir build
cd build
#cmake -DCMAKE_CXX_FLAGS="-std=c++17 " -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$(pwd)/../../PREREQS_INSTALL ..
cmake -DCMAKE_CXX_FLAGS="-std=c++17 " -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(pwd)/../../PREREQS_INSTALL ..
make
make install 
cd ../..

echo
echo "Installing IFOPT."
echo

export IPOPT_DIR=/share/apps/ipopt/3.14/
export EIGEN3_ROOT=$(pwd)/PREREQS_INSTALL/include/eigen3
export EIGEN3_ROOT_DIR=$(pwd)/PREREQS_INSTALL/include/eigen3

cd ifopt 
mkdir build
cd build
#cmake -DCMAKE_CXX_FLAGS="-std=c++17 " -DCMAKE_INSTALL_PREFIX=$(pwd)/../../PREREQS_INSTALL ..
cmake -DCMAKE_CXX_FLAGS="-std=c++17 " -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(pwd)/../../PREREQS_INSTALL ..
make
make install
cd ../..

