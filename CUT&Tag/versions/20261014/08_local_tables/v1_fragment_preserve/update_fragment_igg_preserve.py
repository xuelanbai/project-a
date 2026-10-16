from copy import copy
import math
import os
from pathlib import Path
import re
import tempfile
from zipfile import ZipFile

from openpyxl import load_workbook
from openpyxl.styles import Alignment


base = Path(r"C:\Users\Luca\Documents\Codex\2026-09-17\cut-tag\outputs\01a0d290-1cba-7113-887a-16c91462a481")
source = base / "工作簿1_已填CUTTag统计.xlsx"
output = base / "工作簿1_补充片段长度与IgG对照.xlsx"
fragment_file = Path(r"C:\Users\Luca\Documents\Codex\2026-09-17\cut-tag\fragment_length.tsv")

lines = fragment_file.read_text(encoding="utf-8-sig").splitlines()
if not lines or lines[0] != "sample\tmean_insert_size_bp":
    raise ValueError("Unexpected fragment length TSV header")
lengths = {}
for line in lines[1:]:
    parts = line.split("\t")
    if len(parts) != 2:
        raise ValueError(f"Malformed fragment length row: {line!r}")
    sample_id, raw_value = parts
    value = float(raw_value)
    if sample_id in lengths or not math.isfinite(value) or value <= 0:
        raise ValueError(f"Duplicate sample or invalid length: {line!r}")
    lengths[sample_id] = value

workbook = load_workbook(source)
sheet = workbook["Sheet3"]
sample_names = [sheet[f"C{row}"].value for row in range(3, 13)]
expected = {f"3D4-{group}-{kind}" for group in (119, 137) for kind in
            ("IgG", "Pol2-rep1", "Pol2-rep2", "Tag-rep1", "Tag-rep2")}
assert set(sample_names) == expected and len(sample_names) == 10
expected_ids = {sample.replace("-", "_").replace("_rep1", "_R1").replace("_rep2", "_R2")
                + ("_R1" if sample.endswith("-IgG") else "") for sample in sample_names}
if set(lengths) != expected_ids:
    raise ValueError(f"Fragment length sample mismatch: missing={expected_ids - set(lengths)}, "
                     f"unexpected={set(lengths) - expected_ids}")

sheet["I2"] = "平均插入片段长度（bp）"
sheet["I2"].alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
sheet.column_dimensions["I"].width = 21

for col, title, width in (("AO", "IgG对照样本（按组匹配）", 27),
                          ("AP", "建议分析名称（MACS2）", 48)):
    sheet.column_dimensions[col].width = width
    for row in range(2, 16):
        template = sheet[f"AN{row}"]
        target = sheet[f"{col}{row}"]
        target._style = copy(template._style)
        target.protection = copy(template.protection)
        target.alignment = Alignment(horizontal="center" if row == 2 else "left",
                                     vertical="center", wrap_text=row == 2)
    sheet[f"{col}2"] = title
sheet.column_dimensions["AN"].width = 57

for row, sample in enumerate(sample_names, start=3):
    sample_id = sample.replace("-", "_")
    if sample.endswith("-IgG"):
        sample_id += "_R1"
    else:
        sample_id = sample_id.replace("_rep1", "_R1").replace("_rep2", "_R2")
    sheet[f"I{row}"] = lengths[sample_id]
    sheet[f"I{row}"].number_format = "0.0"
    group = sample.split("-")[1]
    igg = f"3D4-{group}-IgG"
    if sample == igg:
        sheet[f"AO{row}"] = "不适用（IgG对照）"
        sheet[f"AP{row}"] = f"{sample}_IgG-control"
    else:
        sheet[f"AO{row}"] = igg
        sheet[f"AP{row}"] = f"{sample}_vs_{igg}_MACS2"

workbook.save(output)

# openpyxl preserves formulas but clears their cached displayed values. Restore
# the unchanged formula caches from the supplied workbook without editing formulas.
cell_pattern = re.compile(rb'<c\b[^>]*\br="([A-Z]+[0-9]+)"[^>]*(?<!/)>.*?</c>', re.S)
value_pattern = re.compile(rb'<v>.*?</v>', re.S)
with ZipFile(source) as original_zip:
    source_caches = {}
    for name in original_zip.namelist():
        if not name.startswith('xl/worksheets/') or not name.endswith('.xml'):
            continue
        caches = {}
        for match in cell_pattern.finditer(original_zip.read(name)):
            value = value_pattern.search(match.group())
            if b'<f' in match.group() and value:
                caches[match.group(1)] = value.group()
        source_caches[name] = caches
print('SOURCE_CACHES', {name: len(caches) for name, caches in source_caches.items()})

handle = tempfile.NamedTemporaryFile(prefix='fragment_igg_', suffix='.xlsx',
                                     dir=output.parent, delete=False)
temporary = Path(handle.name)
handle.close()
restored = 0
restored_keys = set()
try:
    with ZipFile(output) as generated_zip, ZipFile(temporary, 'w') as fixed_zip:
        for info in generated_zip.infolist():
            content = generated_zip.read(info.filename)
            caches = source_caches.get(info.filename)
            if caches:
                def restore(match):
                    global restored
                    cell = match.group()
                    cached = caches.get(match.group(1))
                    if cached is None or b'<f' not in cell:
                        return cell
                    restored += 1
                    restored_keys.add((info.filename, match.group(1)))
                    if value_pattern.search(cell):
                        return value_pattern.sub(cached, cell, count=1)
                    return cell.replace(b'</c>', cached + b'</c>')
                content = cell_pattern.sub(restore, content)
            fixed_zip.writestr(info, content)
    if restored != sum(len(caches) for caches in source_caches.values()):
        missing = [(name, addr.decode()) for name, caches in source_caches.items()
                   for addr in caches if (name, addr) not in restored_keys]
        print('MISSING_CACHES', missing)
        raise RuntimeError(f'Only restored {restored} cached formula values')
    os.replace(temporary, output)
finally:
    if temporary.exists():
        temporary.unlink()

print(f'RESTORED_FORMULA_CACHES={restored}')
print(output.as_posix())
