#!/bin/bash

echo
echo "Building metabolism"
echo

cd $PROJECT_DIR

source bash_scripts/env.sh

mkdir build
cd build
cmake ..
make