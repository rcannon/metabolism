#!/bin/sh
module purge
module load gcc/9.1.0
module load cmake/3.21.4

# help IFOPT find IPOPT
export IPOPT_DIR=$(pwd)/Ipopt
export EIGEN3_ROOT=$(pwd)/include/eigen3
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$IPOPT_DIR/lib

#TODO: wait for stack exchange up to add files to #include search path