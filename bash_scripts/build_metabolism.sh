#!/bin/bash

echo
echo "Building metabolism"
echo

PROJECT_DIR=$(dirname $0)/..
cd $PROJECT_DIR

mkdir build
cd build
cmake ..
make