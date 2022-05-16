#!/bin/bash

echo
echo "Building metabolism"
echo

PROJECT_DIR=$(dirname $0)/..
cd $PROJECT_DIR

source bash_scripts/env.sh

mkdir build
cd build
cmake ..
make