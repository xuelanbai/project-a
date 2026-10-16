#!/usr/bin/env bash
set -euo pipefail
bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
analysis_env="/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus"
bedtools_env="/mnt/nas/swine/home/jxwu/Projects/FOODIE/envs/bedtools_intersect"
for path in "$bundle_dir/config.json" "$bundle_dir/run.py" "$analysis_env/bin/python" \
            "$bedtools_env/bin/bedtools"; do
    [[ -s "$path" ]] || { echo "MISSING_OR_EMPTY $path" >&2; exit 1; }
done
"$analysis_env/bin/python" - <<'PY'
import importlib.util
for module in ("pysam", "matplotlib"):
    if importlib.util.find_spec(module) is None:
        raise SystemExit(f"missing Python module: {module}")
print("python_modules=OK")
PY
"$bedtools_env/bin/bedtools" --version
echo "environment=OK"
