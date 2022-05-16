#!/bin/sh

echo
echo "cleaning..."
echo

PROJECT_DIR=$(dirname $0)/..
cd $PROJECT_DIR

rm -rf build
rm -rf eigen/build
rm -rf ifopt/build
rm -rf PREREQS_INSTALL