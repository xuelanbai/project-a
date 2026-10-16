#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
config="${1:-${script_dir}/config.json}"
python_bin="/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus/bin/python"
bedtools_bin="/mnt/nas/swine/home/jxwu/Projects/FOODIE/envs/bedtools_intersect/bin/bedtools"

[[ -x "${python_bin}" ]] || { echo "Missing Python: ${python_bin}" >&2; exit 1; }
[[ -x "${bedtools_bin}" ]] || { echo "Missing BEDTools: ${bedtools_bin}" >&2; exit 1; }
[[ -s "${config}" ]] || { echo "Missing config: ${config}" >&2; exit 1; }
"${python_bin}" -c 'import pysam; print("pysam", pysam.__version__)'
"${bedtools_bin}" --version
"${python_bin}" -m py_compile "${script_dir}/run.py" "${script_dir}/test_small.py"
"${python_bin}" "${script_dir}/test_small.py"
echo "ENVIRONMENT_AND_SMALL_TEST_OK"
