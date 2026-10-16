#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
config="${1:-${script_dir}/config.json}"
if [[ "${config}" != /* ]]; then
    echo 'Config path must be absolute (or omit it to use the bundled config).' >&2
    exit 1
fi
cd "${script_dir}"
exec sbatch --export="ALL,CONFIG=${config}" "${script_dir}/run.slurm"
