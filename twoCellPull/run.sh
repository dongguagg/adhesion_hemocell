#!/usr/bin/env bash
set -euo pipefail

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
set +e
source "${case_dir}/../src_hemocell/loadHemoCell.sh"
environment_status=$?
set -e

if [ "${environment_status}" -ne 0 ]; then
    echo "WARNING: optional HemoCell environment checks failed; continuing because mpirun is available." >&2
fi

if ! command -v mpirun >/dev/null 2>&1; then
    echo "ERROR: mpirun was not found" >&2
    exit 1
fi

cd "${case_dir}"
mpirun -np 12 ./twoCellPull config.xml 2>&1 | tee out.txt
