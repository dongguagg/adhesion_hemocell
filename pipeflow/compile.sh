#!/usr/bin/env bash
set -eo pipefail

case_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
hemocell_source_dir="$(cd -- "$case_dir/../src_hemocell" && pwd)"
source "$hemocell_source_dir/loadHemoCell.sh"

cmake -S "$hemocell_source_dir" -B "$hemocell_source_dir/build" \
    -DBUILD_TESTING=OFF
cmake --build "$hemocell_source_dir/build" --target hemocell \
    --parallel "${JOBS:-2}"
test -f "$hemocell_source_dir/build/libhemocell.a"

cmake -S "$case_dir" -B "$case_dir/build" \
    -DHEMOCELL_SOURCE_DIR="$hemocell_source_dir"
cmake --build "$case_dir/build" --target pipeflow --parallel "${JOBS:-2}"
test -x "$case_dir/pipeflow"
