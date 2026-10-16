"""Annotate host candidate genes and test full TF gene sets against callable genes.

The universe is supplied independently for each state and must contain every
candidate gene in that state. Gene IDs match exactly; this module never guesses
symbols, strips ID versions, or expands orthologs. GO BP and Reactome each use
their own annotated subset of that state's callable-gene universe.

Enrichment uses every candidate gene set before theme-based shortlisting. All
background-size-eligible terms, including zero-hit terms, enter the testing
family. BH adjustment is calculated both within TF/state/source and across all
TF/state/source/term tests within each group_level. Motif-name and program
analyses remain separate. An association or a theme match does not establish TF
activity, direct regulation, causality, or a difference between states. The
shortlist is explicitly a functional-theme selection, not a confidence ranking.

The annotation TSV must already contain primary Sus scrofa gene identifiers.
Evidence, URLs and annotation releases are retained without external lookups.
SciPy is required only when a nonempty term family is tested.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path


ANNOTATION_COLUMNS = (
    "gene_id", "source", "term_id", "term_name", "category", "evidence",
    "source_url", "annotation_release", "species",
)
CANDIDATE_COLUMNS = ("contrast", "gene_id", "gene_symbol", "group_level", "tf_label")
SOURCES = ("GO_BP", "Reactome")
CATEGORIES = {"lipid_metabolism", "inflammation_immune", "both", "other"}
GROUP_LEVELS = {"motif_name", "program"}
FUNCTION_COLUMNS = (
    "candidate_row_index", "function_annotation_status", "function_sources",
    "function_primary_category", "function_categories", "function_term_count", "function_terms_json",
    "function_lipid_metabolism", "function_inflammation_immune",
    "lipid_term_ids", "inflammation_term_ids",
    "function_selection_basis",
)
EVIDENCE_COLUMNS = (
    "candidate_row_index", "contrast", "group_level", "tf_label", "gene_symbol",
) + ANNOTATION_COLUMNS
ENRICHMENT_COLUMNS = (
    "contrast", "group_level", "tf_label", "source", "term_id", "term_name",
    "category", "k", "n", "K", "N", "enrichment_fold", "odds_ratio",
    "p_value", "p_adjust_bh_tf_state_source", "p_adjust_bh_global_group_level",
    "hit_gene_ids", "annotation_releases",
)
COVERAGE_COLUMNS = (
    "contrast", "group_level", "tf_label", "source", "raw_background_genes",
    "raw_candidate_genes", "annotated_background_genes",
    "unannotated_background_excluded", "annotated_candidate_genes",
    "unannotated_candidate_excluded", "background_terms",
    "terms_below_min_size", "terms_above_max_size", "tested_terms",
    "zero_hit_terms", "status",
)


@contextmanager
def _open_text(path, mode="rt"):
    """Read ordinary or gzipped TSVs without changing their identifiers."""
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, mode, encoding="utf-8-sig" if "r" in mode else "utf-8",
                newline="") as handle:
        yield handle


def load_annotations(annotation_path):
    """Validate the normalized annotation contract; collapse identical rows.

    Multiple evidence records may support one gene-term membership. They remain
    available in the evidence output but never inflate gene or term counts.
    ``None`` explicitly means that annotations have not been provided.
    """
    if annotation_path is None:
        return [], {}
    records, seen, term_metadata = [], set(), {}
    with _open_text(annotation_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        header = reader.fieldnames
        if (header is None or len(header) != len(set(header))
                or set(header) != set(ANNOTATION_COLUMNS)):
            raise ValueError("Annotation header must contain exactly these unique columns: "
                             + ", ".join(ANNOTATION_COLUMNS))
        try:
            for line_number, raw in enumerate(reader, start=2):
                if None in raw or any(value is None for value in raw.values()):
                    raise ValueError(f"Malformed annotation row at line {line_number}")
                row = {key: raw[key].strip() for key in ANNOTATION_COLUMNS}
                if any(not value for value in row.values()):
                    raise ValueError(f"Blank annotation value at line {line_number}")
                if any("\n" in value or "\r" in value or "\t" in value
                       for value in row.values()):
                    raise ValueError(f"Embedded control character at line {line_number}")
                if row["source"] not in SOURCES:
                    raise ValueError(f"Unsupported annotation source at line {line_number}: "
                                     + row["source"])
                if row["category"] not in CATEGORIES:
                    raise ValueError(f"Unsupported annotation category at line {line_number}: "
                                     + row["category"])
                if row["species"] != "Sus scrofa":
                    raise ValueError(f"Expected Sus scrofa annotation at line {line_number}")
                key = (row["source"], row["term_id"])
                metadata = (row["term_name"], row["category"])
                if key in term_metadata and term_metadata[key] != metadata:
                    raise ValueError(f"Conflicting term name/category for {key!r}")
                term_metadata[key] = metadata
                identity = tuple(row[column] for column in ANNOTATION_COLUMNS)
                if identity not in seen:
                    records.append(row)
                    seen.add(identity)
        except csv.Error as exc:
            raise ValueError(f"Malformed annotation TSV: {exc}") from exc
    return records, term_metadata


def _bh(p_values):
    """Benjamini-Hochberg adjusted p-values in input order, including p=1."""
    size = len(p_values)
    adjusted = [1.0] * size
    ceiling = 1.0
    for rank, index in reversed(list(enumerate(sorted(range(size),
                                                    key=p_values.__getitem__), 1))):
        ceiling = min(ceiling, p_values[index] * size / rank)
        adjusted[index] = ceiling
    return adjusted


def _hypergeom_sf(k, n, K, N):
    try:
        from scipy.stats import hypergeom
    except ImportError as exc:
        raise RuntimeError("Pathway enrichment requires scipy in the server environment") from exc
    return float(hypergeom.sf(k - 1, N, K, n))


def _odds_ratio(k, n, K, N):
    """Foreground versus the disjoint universe-minus-foreground contingency."""
    numerator = k * (N - K - n + k)
    denominator = (n - k) * (K - k)
    if denominator:
        return numerator / denominator
    return math.inf if numerator else None


def _write_tsv(path, rows, columns):
    # Exclusive creation also protects against a file appearing after preflight.
    with _open_text(path, "xt") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: "NA" if row.get(column) is None else row.get(column, "")
                             for column in columns})


def _validate_candidates(candidate_rows, universes):
    rows, candidate_columns = [], list(CANDIDATE_COLUMNS)
    raw_universes = {}
    for state, genes in universes.items():
        if not isinstance(state, str) or not state.strip() or state != state.strip():
            raise ValueError("Universe state keys must be nonempty, unpadded strings")
        if isinstance(genes, str):
            raise ValueError("Each universe must be a collection of gene IDs, not a string")
        raw_universes[state] = set(genes)
        if any(not isinstance(gene, str) or not gene.strip() or gene != gene.strip()
               for gene in raw_universes[state]):
            raise ValueError(f"Universe {state!r} contains a blank or malformed gene ID")
    for index, original in enumerate(candidate_rows, 1):
        row = dict(original)
        missing = set(CANDIDATE_COLUMNS) - set(row)
        if missing:
            raise ValueError(f"Candidate row {index} missing fields: {sorted(missing)}")
        if set(row) & set(FUNCTION_COLUMNS):
            raise ValueError(f"Candidate row {index} uses reserved function-output fields")
        for field in ("contrast", "gene_id", "group_level", "tf_label"):
            value = row[field]
            if not isinstance(value, str) or not value.strip() or value != value.strip():
                raise ValueError(f"Candidate row {index} has invalid {field}")
        if row["group_level"] not in GROUP_LEVELS:
            raise ValueError(f"Candidate row {index} has unsupported group_level")
        state, gene = row["contrast"], row["gene_id"]
        if state not in raw_universes:
            raise ValueError(f"No callable-gene universe for candidate state {state!r}")
        if gene not in raw_universes[state]:
            raise ValueError(f"Candidate gene {gene!r} is outside raw universe {state!r}")
        for column in row:
            if column not in candidate_columns:
                candidate_columns.append(column)
        rows.append(row)
    return rows, raw_universes, candidate_columns


def run_enrichment(candidate_rows, universes, annotation_path, outdir,
                   min_term_size=5, max_term_size=1000):
    """Write row-preserving annotations, evidence, full enrichment and shortlist.

    ``candidate_rows`` is an iterable of dictionaries with contrast, gene_id,
    gene_symbol, group_level (motif_name/program), and tf_label. Additional fields
    are preserved. ``universes`` maps exactly those contrast labels to all
    callable linked gene IDs, independently per state. ``annotation_path`` is a
    normalized TSV/TSV.gz or None. Bounds are inclusive and apply to the number
    of distinct annotated callable genes in a term, separately in each state.

    A zero-annotated-foreground group retains eligible hypotheses with p=1 and
    undefined fold/odds; a missing annotated background produces coverage only.
    All outputs include headers even when empty. Invalid inputs and existing
    output paths fail before any output files are written; files are also opened
    with exclusive creation. Returns JSON-serializable summary counts and paths.
    """
    if (isinstance(min_term_size, bool) or isinstance(max_term_size, bool)
            or not isinstance(min_term_size, int) or not isinstance(max_term_size, int)
            or min_term_size < 1 or max_term_size < min_term_size):
        raise ValueError("Term-size bounds must be integers with 1 <= min <= max")
    outdir = Path(outdir)
    paths = {
        "candidate_functions": outdir / "04_TF_candidate_functions.tsv",
        "gene_function_evidence": outdir / "04_gene_function_evidence.tsv.gz",
        "pathway_enrichment": outdir / "05_TF_pathway_enrichment.tsv",
        "enrichment_coverage": outdir / "05_enrichment_coverage.tsv",
        "candidate_shortlist": outdir / "06_candidate_shortlist.tsv",
    }
    existing = [str(path) for path in paths.values() if path.exists() or path.is_symlink()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing enrichment output: " + ", ".join(existing))
    if outdir.exists() and not outdir.is_dir():
        raise FileExistsError(f"Enrichment output directory is an existing non-directory: {outdir}")
    candidates, raw_universes, candidate_columns = _validate_candidates(candidate_rows, universes)
    annotations, term_metadata = load_annotations(annotation_path)
    by_gene, source_genes, term_genes, term_releases = (
        defaultdict(list), defaultdict(set), defaultdict(set), defaultdict(set)
    )
    for annotation in annotations:
        gene, source, term = annotation["gene_id"], annotation["source"], annotation["term_id"]
        by_gene[gene].append(annotation)
        source_genes[source].add(gene)
        term_genes[source, term].add(gene)
        term_releases[source, term].add(annotation["annotation_release"])
    for records in by_gene.values():
        records.sort(key=lambda row: tuple(row[column] for column in ANNOTATION_COLUMNS))

    functions, shortlist, groups = [], [], defaultdict(set)
    for index, row in enumerate(candidates, 1):
        annotations_for_gene = by_gene.get(row["gene_id"], [])
        unique_terms = sorted({(a["source"], a["term_id"], a["term_name"], a["category"])
                               for a in annotations_for_gene})
        categories = {term[3] for term in unique_terms}
        lipid = bool(categories & {"lipid_metabolism", "both"})
        immune = bool(categories & {"inflammation_immune", "both"})
        primary_category = ("both" if lipid and immune else "lipid_metabolism" if lipid
                            else "inflammation_immune" if immune else "other" if unique_terms
                            else "unannotated")
        enriched_row = dict(row)
        enriched_row.update({
            "candidate_row_index": index,
            "function_annotation_status": "annotated" if unique_terms else "unannotated",
            "function_sources": ";".join(sorted({term[0] for term in unique_terms})),
            "function_primary_category": primary_category,
            "function_categories": ";".join(sorted(categories)),
            "function_term_count": len(unique_terms),
            "function_terms_json": json.dumps([
                dict(zip(("source", "term_id", "term_name", "category"), term))
                for term in unique_terms], ensure_ascii=False, separators=(",", ":")),
            "function_lipid_metabolism": int(lipid),
            "function_inflammation_immune": int(immune),
            "lipid_term_ids": ";".join(f"{source}:{term}" for source, term, _, category in unique_terms
                                       if category in {"lipid_metabolism", "both"}),
            "inflammation_term_ids": ";".join(f"{source}:{term}" for source, term, _, category in unique_terms
                                              if category in {"inflammation_immune", "both"}),
            "function_selection_basis": (
                "functional_theme_annotation_match_only;not_regulatory_confidence"
                if lipid or immune else ""),
        })
        functions.append(enriched_row)
        if lipid or immune:
            shortlist.append(enriched_row)
        groups[row["contrast"], row["group_level"], row["tf_label"]].add(row["gene_id"])

    enrichment, coverage = [], []
    background_cache = {}
    for state, raw in sorted(raw_universes.items()):
        for source in SOURCES:
            annotated = raw & source_genes[source]
            terms = {term: genes & annotated for (term_source, term), genes in term_genes.items()
                     if term_source == source and genes & annotated}
            eligible = {term: genes for term, genes in terms.items()
                        if min_term_size <= len(genes) <= max_term_size}
            background_cache[state, source] = annotated, terms, eligible

    for (state, group_level, tf_label), raw_foreground in sorted(groups.items()):
        for source in SOURCES:
            background, terms, eligible = background_cache[state, source]
            foreground = raw_foreground & background
            N, n = len(background), len(foreground)
            zero_hits = 0
            for term, members in sorted(eligible.items()):
                hits = foreground & members
                k, K = len(hits), len(members)
                zero_hits += int(k == 0)
                name, category = term_metadata[source, term]
                enrichment.append({
                    "contrast": state, "group_level": group_level, "tf_label": tf_label,
                    "source": source, "term_id": term, "term_name": name,
                    "category": category, "k": k, "n": n, "K": K, "N": N,
                    "enrichment_fold": (k / n) / (K / N) if n else None,
                    "odds_ratio": _odds_ratio(k, n, K, N),
                    "p_value": _hypergeom_sf(k, n, K, N),
                    "hit_gene_ids": ";".join(sorted(hits)),
                    "annotation_releases": ";".join(sorted(term_releases[source, term])),
                })
            if not N:
                status = "no_annotated_background"
            elif not n:
                status = "no_annotated_foreground"
            elif not eligible:
                status = "no_terms_in_size_range"
            else:
                status = "tested"
            coverage.append({
                "contrast": state, "group_level": group_level, "tf_label": tf_label,
                "source": source, "raw_background_genes": len(raw_universes[state]),
                "raw_candidate_genes": len(raw_foreground),
                "annotated_background_genes": N,
                "unannotated_background_excluded": len(raw_universes[state]) - N,
                "annotated_candidate_genes": n,
                "unannotated_candidate_excluded": len(raw_foreground) - n,
                "background_terms": len(terms),
                "terms_below_min_size": sum(len(genes) < min_term_size for genes in terms.values()),
                "terms_above_max_size": sum(len(genes) > max_term_size for genes in terms.values()),
                "tested_terms": len(eligible), "zero_hit_terms": zero_hits, "status": status,
            })

    within_families, global_families = defaultdict(list), defaultdict(list)
    for index, row in enumerate(enrichment):
        within_families[row["group_level"], row["contrast"], row["tf_label"], row["source"]].append(index)
        global_families[row["group_level"]].append(index)
    for families, field in (
        (within_families, "p_adjust_bh_tf_state_source"),
        (global_families, "p_adjust_bh_global_group_level"),
    ):
        for indexes in families.values():
            adjusted = _bh([enrichment[index]["p_value"] for index in indexes])
            for index, value in zip(indexes, adjusted):
                enrichment[index][field] = value

    def evidence_rows():
        for index, candidate in enumerate(candidates, 1):
            for annotation in by_gene.get(candidate["gene_id"], []):
                yield dict(annotation, candidate_row_index=index,
                           **{column: candidate[column] for column in
                              ("contrast", "group_level", "tf_label", "gene_symbol")})

    outdir.mkdir(parents=True, exist_ok=True)
    output_candidate_columns = candidate_columns + list(FUNCTION_COLUMNS)
    _write_tsv(paths["candidate_functions"], functions, output_candidate_columns)
    _write_tsv(paths["gene_function_evidence"], evidence_rows(), EVIDENCE_COLUMNS)
    _write_tsv(paths["pathway_enrichment"], enrichment, ENRICHMENT_COLUMNS)
    _write_tsv(paths["enrichment_coverage"], coverage, COVERAGE_COLUMNS)
    _write_tsv(paths["candidate_shortlist"], shortlist, output_candidate_columns)
    return {
        "candidate_rows": len(candidates),
        "candidate_groups": len(groups),
        "annotated_candidate_rows": sum(row["function_annotation_status"] == "annotated"
                                        for row in functions),
        "unannotated_candidate_rows": sum(row["function_annotation_status"] == "unannotated"
                                          for row in functions),
        "annotation_evidence_records": len(annotations),
        "candidate_evidence_rows": sum(len(by_gene.get(row["gene_id"], [])) for row in candidates),
        "shortlist_rows": len(shortlist),
        "tested_hypotheses": len(enrichment),
        "tested_hypotheses_by_group_level": {level: len(global_families.get(level, []))
                                             for level in sorted(GROUP_LEVELS)},
        "coverage_rows": len(coverage),
        "annotation_provided": annotation_path is not None,
        "available_sources": sorted({row["source"] for row in annotations}),
        "outputs": {key: str(path) for key, path in paths.items()},
    }
