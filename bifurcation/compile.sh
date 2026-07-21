#!/bin/bash

source ../src_hemocell/loadHemoCell.sh
mkdir -p build
cd build || exit 1
cmake ..
cmake --build .
