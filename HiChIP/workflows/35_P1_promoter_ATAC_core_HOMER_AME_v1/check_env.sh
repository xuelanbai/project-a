#!/usr/bin/env bash
set -euo pipefail
bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
config="${1:-$bundle_dir/config.json}"
analysis_env="/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus"
motif_env="/mnt/nas/swine/home/jxwu/Projects/A/envs/motif"

for path in "$config" "$bundle_dir/run.py" "$analysis_env/bin/python" \
            "$analysis_env/bin/findMotifs.pl" "$motif_env/bin/ame"; do
    [[ -s "$path" ]] || { echo "MISSING_OR_EMPTY $path" >&2; exit 1; }
done
"$analysis_env/bin/python" - <<'PY'
import importlib.util
for module in ("pysam", "matplotlib"):
    if importlib.util.find_spec(module) is None:
        raise SystemExit(f"missing Python module: {module}")
print("python_modules=OK")
PY
"$motif_env/bin/ame" --version | head -n 1
"$analysis_env/bin/findMotifs.pl" 2>&1 | head -n 1 || true
echo "environment=OK"
