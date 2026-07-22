#!/bin/bash

source ../../loadHemoCell.sh
set -e

cmake -S . -B build
cmake --build build --parallel "$(nproc)"
