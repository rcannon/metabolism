#!/bin/bash

echo
echo "Building metabolism"
echo

cd $PROJECT_DIR

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make