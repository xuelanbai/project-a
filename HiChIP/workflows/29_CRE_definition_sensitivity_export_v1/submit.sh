#!/usr/bin/env bash
set -euo pipefail

bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
config="${1:-$bundle_dir/config.json}"
if [[ "$bundle_dir" != /* ]]; then
    echo "Bundle path must be absolute." >&2
    exit 1
fi
if [[ "$config" != /* ]]; then
    echo "Config path must be absolute (or omit it)." >&2
    exit 1
fi
exec sbatch --chdir="$bundle_dir" --export="ALL,BUNDLE_DIR=$bundle_dir,CONFIG=$config" "$bundle_dir/run.slurm"
