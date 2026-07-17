#!/usr/bin/env bash

set -euo pipefail
trap 'exit 130' INT

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
HEMOCELL_SOURCE_DIR="${PROJECT_ROOT}/src_hemocell"
HEMOCELL_BUILD_DIR="${HEMOCELL_SOURCE_DIR}/build"
CASE_BUILD_DIR="${SCRIPT_DIR}/build"
HEMOCELL_ENV="${HEMOCELL_SOURCE_DIR}/loadHemoCell.sh"
HEMOCELL_LIBRARY="${HEMOCELL_BUILD_DIR}/libhemocell.a"
CASE_EXECUTABLE="${SCRIPT_DIR}/twoCellShear"

echo "=========== Building twoCellShear ==========="
date

if [[ ! -f "${HEMOCELL_ENV}" ]]; then
    echo "ERROR: HemoCell environment script not found: ${HEMOCELL_ENV}" >&2
    exit 1
fi

# loadHemoCell.sh also audits Python post-processing packages such as h5py.
# Those packages are not needed for these two C/C++ build stages, so keep the
# loaded compiler/MPI/HDF5 environment and let CMake perform the required checks.
set +e
# shellcheck source=/dev/null
source "${HEMOCELL_ENV}"
environment_status=$?
set -e
if [[ ${environment_status} -ne 0 ]]; then
    echo "WARNING: loadHemoCell.sh reported a missing component; continuing because CMake will validate all build dependencies." >&2
fi

MPI_C_COMPILER="${MPICC:-/bin/mpicc}"
MPI_CXX_COMPILER="${MPICXX:-/bin/mpicxx}"
HEMOCELL_CXX_FLAGS="${CXXFLAGS:--DOMPI_SKIP_MPICXX}"

echo "[1/2] Configuring and building the HemoCell static library"
cmake \
    -S "${HEMOCELL_SOURCE_DIR}" \
    -B "${HEMOCELL_BUILD_DIR}" \
    -DBUILD_TESTING=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${MPI_C_COMPILER}" \
    -DCMAKE_CXX_COMPILER="${MPI_CXX_COMPILER}" \
    -DMPI_C_COMPILER="${MPI_C_COMPILER}" \
    -DMPI_CXX_COMPILER="${MPI_CXX_COMPILER}" \
    -DMPI_CXX_SKIP_MPICXX=TRUE \
    -DCMAKE_CXX_FLAGS="${HEMOCELL_CXX_FLAGS}"

cmake --build "${HEMOCELL_BUILD_DIR}" --target hemocell --parallel

if [[ ! -f "${HEMOCELL_LIBRARY}" ]]; then
    echo "ERROR: HemoCell build did not create ${HEMOCELL_LIBRARY}" >&2
    exit 1
fi

echo "[2/2] Configuring and building the twoCellShear executable"
cmake \
    -S "${SCRIPT_DIR}" \
    -B "${CASE_BUILD_DIR}" \
    -DHEMOCELL_SOURCE_DIR="${HEMOCELL_SOURCE_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${MPI_C_COMPILER}" \
    -DCMAKE_CXX_COMPILER="${MPI_CXX_COMPILER}" \
    -DMPI_C_COMPILER="${MPI_C_COMPILER}" \
    -DMPI_CXX_COMPILER="${MPI_CXX_COMPILER}" \
    -DMPI_CXX_SKIP_MPICXX=TRUE \
    -DCMAKE_CXX_FLAGS="${HEMOCELL_CXX_FLAGS}"

cmake --build "${CASE_BUILD_DIR}" --target twoCellShear --parallel

if [[ ! -x "${CASE_EXECUTABLE}" ]]; then
    echo "ERROR: twoCellShear build did not create ${CASE_EXECUTABLE}" >&2
    exit 1
fi

date
echo "Build complete: ${CASE_EXECUTABLE}"
