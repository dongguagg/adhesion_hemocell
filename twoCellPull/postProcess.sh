#!/usr/bin/env bash

export PATH="/opt/software/anaconda3/bin:$PATH"
unset PYTHONPATH

cd "$(dirname "${BASH_SOURCE[0]}")"
bash ../src_hemocell/scripts/batchPostProcess.sh
