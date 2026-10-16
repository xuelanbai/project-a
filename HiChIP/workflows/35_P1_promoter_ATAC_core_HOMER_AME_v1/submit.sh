#!/usr/bin/env bash
set -euo pipefail
bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
config="${1:-$bundle_dir/config.json}"
if [[ "$config" != /* ]]; then
    echo 'Config path must be absolute (or omit it to use the bundled config).' >&2
    exit 1
fi
exec sbatch --chdir="$bundle_dir" --export="ALL,BUNDLE_DIR=$bundle_dir,CONFIG=$config" "$bundle_dir/run.slurm"
