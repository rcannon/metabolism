#!/bin/sh
module purge
module load gcc/9.1.0
module load cmake/3.21.4

# help IFOPT find IPOPT
#export IPOPT_DIR=
export EIGEN3_ROOT=$(pwd)/include/eigen3

#TODO: wait for stack exchange up to add files to #include search path