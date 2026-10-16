import json
from pathlib import Path

source = Path(r"F:\2026.09.24_CUT&TAG\results\cuttag_table_inputs.txt")
text = source.read_text(encoding="utf-8")
start = text.index("### DATA\n") + len("### DATA\n")
end = text.index("\n## FLAGSTAT", start)
blocks = json.loads(text[start:end])

for i, block in enumerate(blocks):
    print(f"BLOCK {i}")
    for sample, metrics in sorted(block.items()):
        chosen = {
            k: metrics.get(k)
            for k in [
                "insert_size_average",
                "error_rate",
                "raw_total_sequences",
                "reads_mapped",
                "PERCENT_DUPLICATION",
                "overall_alignment_rate",
            ]
            if k in metrics
        }
        if chosen:
            print(sample, chosen)
