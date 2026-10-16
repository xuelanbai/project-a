import json
import re
from pathlib import Path

SOURCE = Path(r"F:\2026.09.24_CUT&TAG\results\cuttag_table_inputs.txt")
text = SOURCE.read_text(encoding="utf-8")

data_start = text.index("### DATA\n") + len("### DATA\n")
data_end = text.index("\n## FLAGSTAT", data_start)
general_blocks = json.loads(text[data_start:data_end])

merged = {}
for block_index, block in enumerate(general_blocks):
    for sample, metrics in block.items():
        rec = merged.setdefault(sample, {})
        for key, value in metrics.items():
            tagged_key = key
            if tagged_key in rec:
                tagged_key = f"block{block_index}_{key}"
            rec[tagged_key] = value

flag_start = text.index("## FLAGSTAT\n") + len("## FLAGSTAT\n")
flag_end = text.index("\n## PICARD_DUPLICATION", flag_start)
flag_text = text[flag_start:flag_end]
for match in re.finditer(r"### (\S+)\.flagstat\n(.*?)(?=\n### |\Z)", flag_text, re.S):
    sample = match.group(1)
    body = match.group(2)
    rec = merged.setdefault(sample, {})
    patterns = {
        "flagstat_total": r"^(\d+) \+ \d+ in total",
        "flagstat_duplicates": r"^(\d+) \+ \d+ duplicates",
        "flagstat_mapped": r"^(\d+) \+ \d+ mapped \(([0-9.]+)%",
        "flagstat_read1": r"^(\d+) \+ \d+ read1",
        "flagstat_read2": r"^(\d+) \+ \d+ read2",
        "flagstat_proper": r"^(\d+) \+ \d+ properly paired \(([0-9.]+)%",
    }
    for key, pattern in patterns.items():
        found = re.search(pattern, body, re.M)
        if found:
            rec[key] = int(found.group(1))
            if found.lastindex and found.lastindex > 1:
                rec[f"{key}_pct"] = float(found.group(2))

picard_start = text.index("## PICARD_DUPLICATION\n") + len("## PICARD_DUPLICATION\n")
picard_end = text.index("\n## PEAK_COUNTS", picard_start)
picard_lines = [line for line in text[picard_start:picard_end].splitlines() if line.strip()]
for line in picard_lines[1:]:
    sample, unique_pairs, optical, nonoptical = line.split("\t")
    rec = merged.setdefault(sample, {})
    rec["picard_unique_pairs"] = float(unique_pairs)
    rec["picard_duplicate_pairs_optical"] = float(optical)
    rec["picard_duplicate_pairs_nonoptical"] = float(nonoptical)
    total_pairs = float(unique_pairs) + float(optical) + float(nonoptical)
    rec["picard_total_pairs"] = total_pairs
    rec["picard_duplicate_pairs"] = float(optical) + float(nonoptical)
    rec["picard_percent_duplication_calc"] = (
        rec["picard_duplicate_pairs"] / total_pairs * 100 if total_pairs else None
    )

peak_start = text.index("## PEAK_COUNTS\n") + len("## PEAK_COUNTS\n")
peak_end = text.index("\n## SEACR_CONSENSUS_COUNTS", peak_start)
peak_lines = [line for line in text[peak_start:peak_end].splitlines() if line.strip()]
for line in peak_lines[1:]:
    sample, caller, peak_type, raw, filtered, path = line.split("\t")
    rec = merged.setdefault(sample, {})
    prefix = caller.lower()
    rec[f"{prefix}_peak_type"] = peak_type
    rec[f"{prefix}_raw_peaks"] = int(raw)
    rec[f"{prefix}_filtered_peaks"] = None if filtered == "n.a." else int(filtered)
    rec[f"{prefix}_peak_path"] = path

cons_start = text.index("## SEACR_CONSENSUS_COUNTS\n") + len("## SEACR_CONSENSUS_COUNTS\n")
cons_end = text.index("\n## TSS_MATRIX_HEADERS", cons_start)
consensus = {}
cons_lines = [line for line in text[cons_start:cons_end].splitlines() if line.strip()]
for line in cons_lines[1:]:
    group, count, path = line.split("\t")
    consensus[group] = {"count": int(count), "path": path}

print("SAMPLE_KEYS")
all_keys = sorted({key for rec in merged.values() for key in rec})
for key in all_keys:
    print(key)

print("\nMERGED_SAMPLE_DATA")
print(json.dumps(dict(sorted(merged.items())), ensure_ascii=False, indent=2))
print("\nCONSENSUS")
print(json.dumps(consensus, ensure_ascii=False, indent=2))
