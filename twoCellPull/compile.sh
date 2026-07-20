#!/usr/bin/env bash
set -euo pipefail

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
hemocell_dir="$(cd "${case_dir}/../src_hemocell" && pwd)"

# loadHemoCell.sh also checks optional Python post-processing packages. The
# benchmark build itself only requires the compiler/MPI/HDF5 toolchain below.
set +e
source "${hemocell_dir}/loadHemoCell.sh"
environment_status=$?
set -e

if [ "${environment_status}" -ne 0 ]; then
    echo "WARNING: optional HemoCell environment checks failed; continuing after checking the C++ build toolchain." >&2
fi

for required in "${MPICC}" "${MPICXX}" cmake; do
    if ! command -v "${required}" >/dev/null 2>&1; then
        echo "ERROR: required build tool not found: ${required}" >&2
        exit 1
    fi
done

cmake -S "${hemocell_dir}" -B "${hemocell_dir}/build" \
    -DBUILD_TESTING=OFF \
    -DCMAKE_C_COMPILER="${MPICC}" \
    -DCMAKE_CXX_COMPILER="${MPICXX}" \
    -DMPI_C_COMPILER="${MPICC}" \
    -DMPI_CXX_COMPILER="${MPICXX}" \
    -DMPI_CXX_SKIP_MPICXX=TRUE \
    -DCMAKE_CXX_FLAGS=-DOMPI_SKIP_MPICXX
cmake --build "${hemocell_dir}/build" --target hemocell --parallel

cmake -S "${case_dir}" -B "${case_dir}/build" \
    -DCMAKE_C_COMPILER="${MPICC}" \
    -DCMAKE_CXX_COMPILER="${MPICXX}" \
    -DMPI_C_COMPILER="${MPICC}" \
    -DMPI_CXX_COMPILER="${MPICXX}" \
    -DMPI_CXX_SKIP_MPICXX=TRUE \
    -DCMAKE_CXX_FLAGS=-DOMPI_SKIP_MPICXX
cmake --build "${case_dir}/build" --target twoCellPull --parallel
