#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NP="${NP:-2}"
CONFIG="${1:-config.xml}"

if (( NP < 2 )); then
    echo "ERROR: the pre-inlet case requires at least two MPI processes." >&2
    exit 1
fi

if ! source "${SCRIPT_DIR}/../src_hemocell/loadHemoCell.sh"; then
    echo "WARNING: optional HemoCell environment checks failed; continuing with the loaded MPI and HDF5 paths." >&2
fi

cd "${SCRIPT_DIR}"
mpirun -np "${NP}" ./bifurcation "${CONFIG}" 2>&1 | tee out.txt
