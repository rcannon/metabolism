#!/bin/sh
module purge
module load gcc/9.1.0
module load cmake/3.21.4

export EIGEN3_ROOT=$(pwd)/include/eigen3
export IPOPT_DIR=$(pwd)/Ipopt

# build METIS
cd ThirdParty-Metis
./get.Metis
./configure --prefix=$IPOPT_DIR
make
make install
cd ..

# build MUMPS
cd ThirdParty-Mumps
./get.Mumps
./configure --with-metis --with-metis-lflags="-L${IPOPT_DIR}/lib -lcoinmetis" \
     --with-metis-cflags="-I${IPOPT_DIR}/include -I${IPOPT_DIR}/include/coin-or -I${IPOPT_DIR}/include/coin-or/metis" \
     --prefix=$IPOPT_DIR CFLAGS="-I${IPOPT_DIR}/include -I${IPOPT_DIR}/include/coin-or -I${IPOPT_DIR}/include/coin-or/metis" \
     FCFLAGS="-I${IPOPT_DIR}/include -I${IPOPT_DIR}/include/coin-or -I${IPOPT_DIR}/include/coin-or/metis"
make
make install
cd ..

# build IPOPT
cd $IPOPT_DIR
mkdir build
cd build
../configure --prefix=${IPOPT_DIR} --disable-java --with-mumps --with-mumps-lflags="-L${IPOPT_DIR}/lib -lcoinmumps" \
     --with-mumps-cflags="-I${IPOPT_DIR}/include/coin-or/mumps"
make
make install