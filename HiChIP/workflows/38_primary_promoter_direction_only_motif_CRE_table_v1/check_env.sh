#!/usr/bin/env bash
set -euo pipefail
bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
config="${1:-$bundle_dir/config.json}"
analysis_env="/mnt/nas/swine/home/jxwu/Projects/A/envs/hichip_locus"
motif_env="/mnt/nas/swine/home/jxwu/Projects/A/envs/motif"

for path in "$config" "$bundle_dir/run.py" "$bundle_dir/motif.py" "$bundle_dir/intervals.py" \
            "$analysis_env/bin/python" "$analysis_env/bin/findMotifs.pl" "$motif_env/bin/ame"; do
    [[ -s "$path" ]] || { echo "MISSING_OR_EMPTY $path" >&2; exit 1; }
done

"$analysis_env/bin/python" - "$config" <<'PY'
import importlib.util, json, pathlib, sys
cfg = json.loads(pathlib.Path(sys.argv[1]).read_text())
paths = [
    pathlib.Path(cfg["source_ATAC_core_run"]) / "02_ATAC_consensus_cores.tsv",
    pathlib.Path(cfg["source_P1_promoter_state_run"]) / "04_P1_promoter_signal_and_state_matrix.tsv",
    pathlib.Path(cfg["source_P1_promoter_state_run"]) / "05_P1_promoter_state_gene_associations.tsv",
    pathlib.Path(cfg["source_mark_matrix_run"]) / "01_enhancer_mark_signal_matrix.tsv",
    pathlib.Path(cfg["source_mark_matrix_run"]) / "02_P1_promoter_mark_signal_matrix.tsv",
    pathlib.Path(cfg["reference_fasta"]), pathlib.Path(cfg["reference_fasta"] + ".fai"),
    pathlib.Path(cfg["motif_database"]), *(pathlib.Path(x) for x in cfg["loops"].values()),
]
bad = [str(p) for p in paths if not p.is_absolute() or not p.is_file() or p.stat().st_size == 0]
if bad:
    raise SystemExit("missing or empty inputs:\n" + "\n".join(bad))
for module in ("pysam", "matplotlib"):
    if importlib.util.find_spec(module) is None:
        raise SystemExit(f"missing Python module: {module}")
print(f"inputs=OK n={len(paths)}")
print("python_modules=OK")
PY
"$motif_env/bin/ame" --version | head -n 1
"$analysis_env/bin/findMotifs.pl" 2>&1 | head -n 1 || true
echo "environment=OK"
