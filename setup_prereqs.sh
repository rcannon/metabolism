#!/bin/sh
module purge
module load gcc/9.1.0
module load cmake/3.21.4

export METAB_DIR=$(pwd)

export EIGEN3_ROOT=$(pwd)/include/eigen3
export IPOPT_DIR=$(pwd)/Ipopt

# build MUMPS
echo ""
echo "Building MUMPS"
echo ""
cd ThirdParty-Mumps
./get.Mumps
./configure --prefix=$IPOPT_DIR CFLAGS="-I${IPOPT_DIR}/include -I${IPOPT_DIR}/include/coin-or" \
     FCFLAGS="-I${IPOPT_DIR}/include -I${IPOPT_DIR}/include/coin-or"
make
make install
cd ..
echo ""
echo "Finished MUMPS"
echo ""

# build IPOPT
echo ""
echo "Building IPOPT"
echo ""
cd $IPOPT_DIR
mkdir build
cd build
../configure --prefix=${IPOPT_DIR} --disable-java --with-mumps --with-mumps-lflags="-L${IPOPT_DIR}/lib -lcoinmumps" \
     --with-mumps-cflags="-I${IPOPT_DIR}/include/coin-or/mumps"
make
make install
cd ../..
echo ""
echo "Finished IPOPT"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$IPOPT_INSTALL_DIR/lib

# for ifopt
export IPOPT_DIR=$(pwd)/Ipopt/build

#build IFOPT
echo ""
echo "Building IFOPT"
echo ""
cd ifopt
mkdir build
cd build
cmake .. -DCMAKE_CXX_FLAGS="-std=c++17 " -DCMAKE_INSTALL_PREFIX=$METAB_DIR
make
make install
echo ""
echo "Finished IFOPT"
echo ""