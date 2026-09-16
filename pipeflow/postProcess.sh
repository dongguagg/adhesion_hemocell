#!/usr/bin/env bash
set -e

export PATH="/opt/software/anaconda3/bin:$PATH"
unset PYTHONPATH

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scripts_dir="$case_dir/../src_hemocell/scripts"

for output_dir in "$case_dir"/output "$case_dir"/output_[0-9]*; do
    [[ -d "$output_dir/hdf5" ]] || continue
    echo "Post-processing $output_dir"
    (
        cd "$output_dir"
        python3 "$scripts_dir/FluidHDF5.py"
        python3 "$scripts_dir/CellHDF5toXMF.py" RBC
    )
done
