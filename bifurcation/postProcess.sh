#!/usr/bin/env bash

export PATH="/opt/software/anaconda3/bin:$PATH"
unset PYTHONPATH

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
postprocess_script="$case_dir/../src_hemocell/scripts/batchPostProcess.sh"

for output_dir in "$case_dir"/output "$case_dir"/output_[0-9]*; do
    [[ -d "$output_dir" ]] || continue
    echo "Post-processing $(basename "$output_dir")"
    (
        cd "$output_dir" || exit 1
        bash "$postprocess_script"
    ) || exit 1
done
