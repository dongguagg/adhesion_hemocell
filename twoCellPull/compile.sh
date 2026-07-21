#!/bin/bash

source ../src_hemocell/loadHemoCell.sh

cmake -S ../src_hemocell -B ../src_hemocell/build -DBUILD_TESTING=OFF
cmake --build ../src_hemocell/build --target hemocell

mkdir -p build
cd build || exit 1
cmake ..
cmake --build .
