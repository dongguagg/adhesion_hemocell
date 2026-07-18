#!/bin/bash

source ./loadHemoCell.sh
mkdir -p build
cd build || exit 1
cmake ..
cmake --build . -j 8
