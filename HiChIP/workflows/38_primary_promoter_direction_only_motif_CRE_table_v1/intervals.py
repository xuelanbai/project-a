"""Reconstruct motif-independent gene universes from original scan windows.

This module performs ordinary interval joins. It is intended to run on the
analysis server; the companion tests use only tiny synthetic inputs.

Coordinates are zero-based, half-open throughout. Gene-feature TSS positions
follow the historical convention: GTF start - 1 on + and GTF end - 1 on -.
Only a literal leading ``chr`` is removed from chromosome names. No gene-ID
version normalization, transcript substitution, peak expansion, time-point
union, or differential-loop inference is performed.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import re
import tempfile
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, TextIO


STATE_GROUPS = {
    "loss": "loss_scan_all_0h_only",
    "gain": "gain_scan_all_6h_only",
}
STATE_TIMES = {"loss": "0h", "gain": "6h"}
REGION_COLUMNS = (
    "region_id", "group", "chrom", "start", "end",
    "passes_concordant_motif_foreground",
)
EVIDENCE_COLUMNS = (
    "state", "region_id", "group", "region_chrom", "region_start",
    "region_end", "gene_id", "gene_symbol", "promoter_chrom",
    "promoter_start", "promoter_end", "gene_tss", "gene_strand",
    "loop_time", "loop_id", "anchor1_chrom", "anchor1_start",
    "anchor1_end", "anchor2_chrom", "anchor2_start", "anchor2_end",
    "loop_source", "loop_line", "contact_count", "loop_p", "loop_q",
    "orientation",
)
_INTEGER = re.compile(r"[0-9]+\Z")
_GTF_ATTRIBUTE = re.compile(r'\s*([^\s;]+)\s+"([^"\r\n]*)"\s*(?:;|$)')


def canonical_chrom(chrom: str) -> str:
    """Strip whitespace and only the literal, case-sensitive ``chr`` prefix."""
    value = str(chrom).strip()
    if value.startswith("chr"):
        value = value[3:]
    if not value or any(char.isspace() for char in value):
        raise ValueError(f"Invalid chromosome name: {chrom!r}")
    return value


def _integer(value: Any, context: str, *, positive: bool = False) -> int:
    text = str(value).strip()
    if not _INTEGER.fullmatch(text):
        raise ValueError(f"{context}: expected a nonnegative integer, got {value!r}")
    result = int(text)
    if positive and result == 0:
        raise ValueError(f"{context}: expected a positive integer")
    return result


def _interval(chrom: str, start: Any, end: Any, context: str) -> tuple[str, int, int]:
    result = (canonical_chrom(chrom), _integer(start, context), _integer(end, context))
    if result[2] <= result[1]:
        raise ValueError(f"{context}: interval end must be greater than start")
    return result


def stable_loop_id(a1chrom: str, start1: int, end1: int,
                   a2chrom: str, start2: int, end2: int) -> str:
    """Return the same coordinate identity for either BEDPE orientation."""
    anchors = sorted((_interval(a1chrom, start1, end1, "loop anchor 1"),
                      _interval(a2chrom, start2, end2, "loop anchor 2")))
    return "|".join(f"{chrom}:{start}-{end}" for chrom, start, end in anchors)


def _open_text(path: str | os.PathLike[str]) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


class RegionMap(dict):
    """A normal region dictionary with optional parsing QC attached."""

    def __init__(self) -> None:
        super().__init__()
        self.qc: dict[str, Any] = {}


def load_regions(path: str | os.PathLike[str], require_concordant: bool = True) -> dict:
    """Read selected scan windows, keyed by ``(loss|gain, region_id)``.

    Original coordinates are retained. Windows must be positive and at most
    500 bp, allowing chromosome-edge clipping. Concordance is the supplied
    binary eligibility flag; the loader never consults motif occurrences.
    """
    regions = RegionMap()
    group_states = {group: state for state, group in STATE_GROUPS.items()}
    counts: Counter = Counter()
    lengths: dict[str, Counter] = {state: Counter() for state in STATE_GROUPS}
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if (not header or len(header) != len(set(header))
                or any(not name for name in header)):
            raise ValueError(f"{path}: empty or duplicate TSV column names")
        missing = set(REGION_COLUMNS) - set(header)
        if missing:
            raise ValueError(f"{path}: missing required region columns: {sorted(missing)}")
        for row in reader:
            context = f"{path}:{reader.line_num}"
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{context}: field count differs from TSV header")
            counts["input_rows"] += 1
            flag = row["passes_concordant_motif_foreground"].strip()
            if flag not in {"0", "1"}:
                raise ValueError(f"{context}: concordant foreground flag must be 0 or 1")
            region_id = row["region_id"].strip()
            if not region_id:
                raise ValueError(f"{context}: empty region_id")
            chrom, start, end = _interval(row["chrom"], row["start"], row["end"], context)
            group = row["group"].strip()
            if group not in group_states:
                counts["excluded_other_group"] += 1
                continue
            if require_concordant and flag != "1":
                counts["excluded_nonconcordant"] += 1
                continue
            if end - start > 500:
                raise ValueError(f"{context}: original scan window exceeds 500 bp")
            state = group_states[group]
            key = (state, region_id)
            parsed = dict(row, state=state, region_id=region_id, group=group,
                          chrom=chrom, start=start, end=end,
                          passes_concordant_motif_foreground=int(flag))
            if key in regions:
                previous = regions[key]
                compared = ("chrom", "start", "end", "group", "passes_concordant_motif_foreground")
                if any(previous[field] != parsed[field] for field in compared):
                    raise ValueError(f"{context}: conflicting duplicate region {key!r}")
                counts["duplicate_identical_rows"] += 1
                continue
            regions[key] = parsed
            lengths[state][end - start] += 1
    regions.qc = dict(counts)
    regions.qc["require_concordant"] = require_concordant
    regions.qc["window_length_counts"] = {
        state: dict(sorted(histogram.items())) for state, histogram in lengths.items()
    }
    regions.qc["selected_regions"] = {
        state: sum(lengths[state].values()) for state in STATE_GROUPS
    }
    return regions


def _load_chrom_sizes(source: Mapping | str | os.PathLike[str]) -> dict[str, int]:
    sizes: dict[str, int] = {}
    originals: dict[str, str] = {}

    def add(chrom: Any, length: Any, context: str) -> None:
        original = str(chrom).strip()
        canonical = canonical_chrom(original)
        if canonical in sizes:
            raise ValueError(
                f"{context}: ambiguous or duplicate chromosome aliases "
                f"{originals[canonical]!r} and {original!r}"
            )
        sizes[canonical] = _integer(length, context, positive=True)
        originals[canonical] = original

    if isinstance(source, Mapping):
        for chrom, length in source.items():
            add(chrom, length, "chromosome sizes")
    else:
        with _open_text(source) as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) != 2:
                    raise ValueError(f"{source}:{number}: expected two chromosome-size fields")
                add(fields[0], fields[1], f"{source}:{number}")
    if not sizes:
        raise ValueError("Chromosome sizes are empty")
    return sizes


@dataclass(frozen=True)
class Promoter:
    chrom: str
    start: int
    end: int
    gene_id: str
    gene_symbol: str
    tss: int
    strand: str


class IntervalIndex:
    """Simple chromosome/bin index; each overlap is returned exactly once."""

    def __init__(self, bin_size: int = 100_000) -> None:
        self.bin_size = bin_size
        self._bins: dict[str, dict[int, list[int]]] = defaultdict(lambda: defaultdict(list))
        self._items: list[tuple[int, int, Any]] = []
        self.qc: dict[str, Any] = {}

    @property
    def chromosomes(self) -> set[str]:
        return set(self._bins)

    def add(self, chrom: str, start: int, end: int, value: Any) -> None:
        item_id = len(self._items)
        self._items.append((start, end, value))
        for bin_id in range(start // self.bin_size, (end - 1) // self.bin_size + 1):
            self._bins[chrom][bin_id].append(item_id)

    def overlaps(self, chrom: str, start: int, end: int) -> Iterable[Any]:
        bins = self._bins.get(canonical_chrom(chrom), {})
        seen: set[int] = set()
        for bin_id in range(start // self.bin_size, (end - 1) // self.bin_size + 1):
            for item_id in bins.get(bin_id, ()):
                if item_id in seen:
                    continue
                seen.add(item_id)
                item_start, item_end, value = self._items[item_id]
                if item_start < end and start < item_end:
                    yield value


def _gtf_attributes(value: str, context: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    position = 0
    while position < len(value):
        if not value[position:].strip():
            break
        match = _GTF_ATTRIBUTE.match(value, position)
        if match is None:
            raise ValueError(f"{context}: malformed GTF attributes near {value[position:position+60]!r}")
        key, item = match.groups()
        item = item.strip()
        if key in {"gene_id", "gene_name"} and key in attributes and attributes[key] != item:
            raise ValueError(f"{context}: conflicting duplicate GTF attribute {key!r}")
        attributes[key] = item
        position = match.end()
    return attributes


def load_promoters(gtf: str | os.PathLike[str], chrom_sizes: Mapping | str | os.PathLike[str],
                   half_width: int = 2000) -> tuple[IntervalIndex, dict[str, str]]:
    """Index GTF gene-feature promoters and return exact gene ID -> symbol."""
    sizes = _load_chrom_sizes(chrom_sizes)
    half_width = _integer(half_width, "promoter half-width", positive=True)
    index = IntervalIndex()
    genes: dict[str, str] = {}
    explicit_symbols: dict[str, str] = {}
    seen_promoters: set[Promoter] = set()
    missing_chroms: Counter = Counter()
    counts: Counter = Counter()
    pending: list[tuple[str, int, int, str, int, str]] = []
    with _open_text(gtf) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            context = f"{gtf}:{number}"
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{context}: expected exactly nine GTF fields")
            start = _integer(fields[3], context, positive=True)
            end = _integer(fields[4], context, positive=True)
            if end < start:
                raise ValueError(f"{context}: GTF end precedes start")
            counts["gtf_feature_rows"] += 1
            if fields[2] != "gene":
                counts["excluded_non_gene_features"] += 1
                continue
            counts["gene_feature_rows"] += 1
            if fields[6] not in {"+", "-"}:
                raise ValueError(f"{context}: gene strand must be + or -")
            attributes = _gtf_attributes(fields[8], context)
            gene_id = attributes.get("gene_id", "")
            if not gene_id:
                raise ValueError(f"{context}: gene feature has no gene_id")
            name = attributes.get("gene_name", "")
            if name:
                if gene_id in explicit_symbols and explicit_symbols[gene_id] != name:
                    raise ValueError(f"{context}: contradictory symbols for gene {gene_id!r}")
                explicit_symbols[gene_id] = name
            genes[gene_id] = explicit_symbols.get(gene_id, gene_id)
            chrom = canonical_chrom(fields[0])
            if chrom not in sizes:
                missing_chroms[chrom] += 1
                counts["gene_features_no_chrom_size"] += 1
                continue
            if end > sizes[chrom]:
                raise ValueError(f"{context}: gene feature exceeds chromosome size")
            tss = start - 1 if fields[6] == "+" else end - 1
            promoter_start = max(0, tss - half_width)
            promoter_end = min(sizes[chrom], tss + half_width + 1)
            pending.append((chrom, promoter_start, promoter_end, gene_id, tss, fields[6]))
    # Resolve names after all gene rows so missing names never mask a later name.
    for chrom, start, end, gene_id, tss, strand in pending:
        promoter = Promoter(chrom, start, end, gene_id, genes[gene_id], tss, strand)
        if promoter in seen_promoters:
            counts["duplicate_identical_promoters"] += 1
            continue
        seen_promoters.add(promoter)
        index.add(chrom, start, end, promoter)
    counts["indexed_promoters"] = len(seen_promoters)
    counts["gene_ids"] = len(genes)
    counts["genes_without_symbol"] = sum(gene_id not in explicit_symbols for gene_id in genes)
    index.qc = dict(counts)
    index.qc["missing_chromosomes"] = dict(sorted(missing_chroms.items()))
    index.promoter_keys = {
        (promoter.gene_id, promoter.chrom, promoter.start, promoter.end, promoter.tss)
        for promoter in seen_promoters
    }
    return index, genes


def _number(value: str, context: str, *, probability: bool = False) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"{context}: invalid numeric value {value!r}") from error
    if not math.isfinite(number) or number < 0 or (probability and number > 1):
        raise ValueError(f"{context}: numeric value outside the permitted range: {value!r}")
    return number


def _loop_paths(value: Any) -> list[Path]:
    if isinstance(value, (str, os.PathLike)):
        return [Path(value)]
    paths = [Path(path) for path in value]
    if not paths:
        raise ValueError("Every state must have at least one loop file")
    return paths


