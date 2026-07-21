#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HEMOCELL_SOURCE_DIR="${SCRIPT_DIR}/../src_hemocell"
BUILD_JOBS="${BUILD_JOBS:-4}"

if ! source "${HEMOCELL_SOURCE_DIR}/loadHemoCell.sh"; then
    echo "WARNING: optional HemoCell environment checks failed; continuing with the loaded compiler and HDF5 paths." >&2
fi

cmake -S "${HEMOCELL_SOURCE_DIR}" \
      -B "${HEMOCELL_SOURCE_DIR}/build" \
      -DBUILD_TESTING=OFF \
      -DMPI_CXX_SKIP_MPICXX=TRUE \
      -DCMAKE_CXX_FLAGS=-DOMPI_SKIP_MPICXX
cmake --build "${HEMOCELL_SOURCE_DIR}/build" \
      --target hemocell \
      --parallel "${BUILD_JOBS}"

cmake -S "${SCRIPT_DIR}" \
      -B "${SCRIPT_DIR}/build" \
      -DMPI_CXX_SKIP_MPICXX=TRUE \
      -DCMAKE_CXX_FLAGS=-DOMPI_SKIP_MPICXX
cmake --build "${SCRIPT_DIR}/build" \
      --target bifurcation \
      --parallel "${BUILD_JOBS}"
