#!/usr/bin/env bash
set -euo pipefail
bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
python_bin=/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus/bin/python
macs3_bin=/mnt/nas/swine/home/jxwu/Projects/FOODIE/envs/macs3_callpeak/bin/macs3
homer_bin=/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus/bin/findMotifs.pl
ame_bin=/mnt/nas/swine/home/jxwu/Projects/A/envs/motif/bin/ame
bedtools_bin=/mnt/nas/swine/home/jxwu/Projects/FOODIE/envs/bedtools_intersect/bin/bedtools

for file in "$python_bin" "$macs3_bin" "$homer_bin" "$ame_bin" "$bedtools_bin"; do
    [[ -x "$file" ]] || { echo "MISSING_EXECUTABLE $file" >&2; exit 1; }
done
[[ -s "$bundle_dir/config.json" ]] || { echo "MISSING_CONFIG $bundle_dir/config.json" >&2; exit 1; }
"$python_bin" - <<'PY'
import importlib
for name in ("pysam", "matplotlib"):
    module = importlib.import_module(name)
    print(name, getattr(module, "__version__", "available"))
PY
"$macs3_bin" --version
"$ame_bin" --version
echo "environment_check=PASS"
