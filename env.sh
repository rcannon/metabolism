#!/bin/sh
module purge
module load gcc/9.1.0
module load cmake/3.21.4

# help IFOPT find IPOPT
export IPOPT_DIR=/share/apps/ipopt/3.14/
export EIGEN3_ROOT=$(pwd)/include/eigen3
export EIGEN3_ROOT_DIR=$(pwd)/include/eigen3
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$IPOPT_DIR/lib

#rm -rf build
#mkdir build
#cd build
#make ..
#TODO: wait for stack exchange up to add files to #include search path