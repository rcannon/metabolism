#!/bin/bash

echo
echo "cleaning..."
echo

cd $PROJECT_DIR

rm -rf build
rm -rf eigen/build
rm -rf ifopt/build
rm -rf PREREQS_INSTALL