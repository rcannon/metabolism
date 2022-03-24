#!/bin/sh
module purge
module load gcc/9.1.0
module load cmake/3.21.4

export METAB_DIR=$(pwd)

export EIGEN3_ROOT=$(pwd)/include/eigen3
export IPOPT_DIR=/share/apps/ipopt/3.14

# for ifopt
export IPOPT_DIR=$(pwd)/Ipopt/build

#build IFOPT
echo ""
echo "Building IFOPT"
echo ""
cd ifopt
mkdir build
cd build
mkdir
ifopt_build
cmake .. -DCMAKE_CXX_FLAGS="-std=c++17 " -DCMAKE_INSTALL_PREFIX=$METAB_DIR/ifopt_build
make
make install
echo ""
echo "Finished IFOPT"
echo ""