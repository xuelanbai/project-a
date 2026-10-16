#!/usr/bin/env bash
set -euo pipefail
bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
resume_from="${1:-}"
config="${2:-$bundle_dir/config.json}"
if [[ -n "$resume_from" && "$resume_from" != /* ]]; then
    echo 'Resume path must be absolute.' >&2
    exit 1
fi
if [[ "$config" != /* ]]; then
    echo 'Config path must be absolute (or omit it to use the bundled config).' >&2
    exit 1
fi
exec sbatch --chdir="$bundle_dir" --export="ALL,BUNDLE_DIR=$bundle_dir,CONFIG=$config,RESUME_FROM=$resume_from" "$bundle_dir/run.slurm"
