#!/usr/bin/env python3
"""Small deterministic tests for coordinate conversion and HiChIP joins."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
import tempfile

import intervals
import run as pipeline
import promoter_definitions


def main() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        fimo = root / "fimo.tsv"
        fimo.write_text(
            "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
            "M1\tTF1\tC1\t51\t58\t+\t10\t1e-6\t0.001\tACGTACGT\n",
            encoding="utf-8",
        )
        region = {
            "region_id": "R1", "chrom": "1", "start": 100, "end": 400,
            "H3K27ac_log2FC": 1.0, "ATAC_log2FC": 0.2, "ATAC_stratum": "ATAC_up",
            "core_ids": ["C1"],
        }
        cores = {"C1": {"core_id": "C1", "chrom": "1", "start": 50, "end": 550,
                         "region_id": "R1", "region": region}}
        motifs = {"M1": {"motif_id": "M1", "motif_name": "TF1", "motif_family": "FAM1",
                          "AME_q": 0.01, "q_field": "AME_BH_q_within_contrast",
                          "enrichment_ratio": 1.5, "fimo_p": 1e-4}}
        sites, site_qc = pipeline.parse_sites(fimo, cores, motifs, root)
        assert site_qc == {"reported_sites": 1, "sites_inside_full_CRE": 1}
        assert sites["R1"]["M1"][0]["start"] == 100
        assert sites["R1"]["M1"][0]["end"] == 108

        promoter_index = intervals.IntervalIndex()
        promoter = intervals.Promoter("1", 900, 1101, "GENE1", "Gene1", 1000, "+")
        promoter_index.add("1", promoter.start, promoter.end, promoter)
        loops = root / "loops.bedpe"
        loops.write_text("chr1\t90\t250\tchr1\t950\t1050\t12\t0.001\t0.005\n", encoding="utf-8")
        destination = root / "state"
        destination.mkdir()
        cfg = {"loops": {"strength": str(loops)}, "loop_q": 0.01,
               "loop_times": {"strength": "6h"}}
        candidates, universe, qc = pipeline.connect_state(
            "strength", {"R1": region}, cores, sites, motifs, promoter_index,
            promoter_definitions.TSS, "not_time_specific", cfg, destination
        )
        assert universe == {"GENE1"}
        assert qc["CRE_loop_gene_links"] == 1
        assert qc["motif_evidence_rows"] == 1
        assert len(candidates) == 1
        assert candidates[0]["unique_CRE_count"] == 1
        assert candidates[0]["strict_site_anchor_unique_CRE_count"] == 1
        with gzip.open(destination / "06_motif_CRE_loop_gene_evidence.tsv.gz", "rt", encoding="utf-8") as handle:
            evidence = list(csv.DictReader(handle, delimiter="\t"))
        assert evidence[0]["any_core_overlaps_CRE_anchor"] == "1"
        assert evidence[0]["any_site_overlaps_CRE_anchor"] == "1"

        selected = pipeline.select_motifs([
            {"contrast": pipeline.CONTRASTS["strength"], "motif_id": "M1", "motif_alt_id": "TF1",
             "motif_family": "FAM1", "enrichment_ratio": "1.5",
             "AME_BH_q_within_contrast": "0.01"},
            {"contrast": pipeline.CONTRASTS["strength"], "motif_id": "M2", "motif_alt_id": "TF2",
             "motif_family": "FAM2", "enrichment_ratio": "1.7",
             "AME_BH_q_within_contrast": "0.20"},
        ], "strength", {"ame_q_field": "AME_BH_q_within_contrast", "ame_q": 0.05})
        assert [row["motif_id"] for row in selected] == ["M1"]

        # P1 passes support anywhere in the full H3 peak. P2 tests support only
        # inside the gene-specific H3-peak/TSS-window intersection.
        gtf = root / "genes.gtf"
        gtf.write_text(
            'chr1\ttest\tgene\t1001\t1200\t.\t+\t.\tgene_id "GENE1"; gene_name "Gene1";\n',
            encoding="utf-8",
        )
        peak_paths = {}
        for time in promoter_definitions.TIMES:
            rows = {
                "H3K27ac": "chr1\t800\t1200\n",
                "ATAC": "chr1\t810\t820\n",
                "H3K4me3": "chr1\t1010\t1020\n",
            }
            for assay, text in rows.items():
                path = root / f"{assay}_{time}.bed"
                path.write_text(text, encoding="utf-8")
                peak_paths[f"{assay}_{time}"] = str(path)
        promoter_out = root / "promoters"
        definitions, definition_qc = promoter_definitions.build_promoter_definitions(
            gtf, {"1": 5000}, 100, peak_paths, promoter_out
        )
        assert list(definitions[promoter_definitions.P1]["0h"].overlaps("1", 800, 1200))
        assert not list(definitions[promoter_definitions.P2]["0h"].overlaps("1", 800, 1200))
        assert (promoter_out / promoter_definitions.P1 / "promoter_consensus_0h_6h.bed").is_file()
        assert definition_qc
    print("SMALL_TESTS_OK")


if __name__ == "__main__":
    main()
