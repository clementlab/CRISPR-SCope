#!/usr/bin/env bash
# Run on an owner or allocated compute node, not a login node.
set -euo pipefail

if [[ "${CONDA_DEFAULT_ENV:-}" != "sc_CRISPRSCope" ]]; then
    printf 'Activate the sc_CRISPRSCope conda environment first.\n' >&2
    exit 2
fi

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log_dir="$HOME/tmp_dir"
scratch_parent="/scratch/general/nfs1/$USER"
mkdir -p "$log_dir" "$scratch_parent"

run_dir="$(mktemp -d "$scratch_parent/crisprscope-integration.XXXXXXXX")"
run_name="$(basename "$run_dir")"
run_link="$log_dir/$run_name"
ln -s "$run_dir" "$run_link"

printf 'Pytest log: %s\n' "$log_dir/$run_name.pytest.log"
printf 'Live pipeline output: %s\n' "$run_link"
printf 'Scratch directory to remove when no longer needed: %s\n' "$run_dir"
printf 'Find the pipeline log with: find -L %q -name settings.txt.log -print\n' "$run_link"

cd "$repo_dir"
python -m pytest -v -m integration \
    --basetemp "$run_dir" \
    tests/test_baf3_minimal_integration.py \
    2>&1 | tee "$log_dir/$run_name.pytest.log"
