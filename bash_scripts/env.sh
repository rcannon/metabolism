#!/bin/sh
module purge
module load gcc/9.1.0
module load cmake/3.21.4

# help IFOPT find IPOPT
#export IPOPT_DIR=/share/apps/ipopt/3.14/
#export EIGEN3_ROOT=$(pwd)/INSTALL/include/eigen3
#export EIGEN3_ROOT_DIR=$(pwd)/INSTALL/include/eigen3


#rm -rf build
#mkdir build
#cd build
#make ..
#TODO: wait for stack exchange up to add files to #include search path