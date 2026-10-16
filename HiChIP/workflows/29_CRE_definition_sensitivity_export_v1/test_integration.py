#!/usr/bin/env python3
"""End-to-end synthetic-data test for BED export, counts and v27 QC."""
from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="CRE_definition_export_") as temporary:
        root = Path(temporary)
        project, source_run, output_root = root / "project", root / "v27", root / "results"
        fai, gtf = root / "genome.fa.fai", root / "genes.gtf.gz"
        write(fai, "".join(f"chr{chrom}\t10000\t0\t0\t0\n" for chrom in [*map(str, range(1, 19)), "X", "Y"]))
        with gzip.open(gtf, "wt", encoding="utf-8") as handle:
            handle.write('chr1\ttest\tgene\t5001\t5500\t.\t+\t.\tgene_id "gene1"; gene_name "GENE1";\n')

        write(source_run / "peaks/ATAC_0h_peaks.narrowPeak", "chr1\t4020\t4060\tA0P\t1\t.\t1\t1\t1\t20\nchr1\t8050\t8120\tA0E\t1\t.\t1\t1\t1\t20\n")
        write(source_run / "peaks/ATAC_6h_peaks.narrowPeak", "chr1\t4070\t4120\tA6P\t1\t.\t1\t1\t1\t20\nchr1\t8150\t8220\tA6E\t1\t.\t1\t1\t1\t20\n")
        write(project / "06.peak/h3ac0.bed", "chr1\t4000\t4100\nchr1\t8000\t8200\n")
        write(project / "06.peak/h3ac6.bed", "chr1\t4050\t4150\nchr1\t8100\t8300\n")
        write(project / "06.peak/h3k40.bed", "chr1\t4040\t4080\n")
        write(project / "06.peak/h3k46.bed", "chr1\t4080\t4130\n")

        write(
            source_run / "03_H3K27ac_consensus_membership.tsv",
            "region_id\tchrom\tsource_time\tpeak_start\tpeak_end\n"
            "H3K27ac_CRE_000001\t1\t0h\t8000\t8200\n"
            "H3K27ac_CRE_000001\t1\t6h\t8100\t8300\n",
        )
        write(
            source_run / "03_CRE_signal_matrix.tsv",
            "chrom\tstart\tend\tregion_id\n1\t8000\t8300\tH3K27ac_CRE_000001\n",
        )

        config = {
            "project": str(project),
            "source_v27_run": str(source_run),
            "output_root": str(output_root),
            "gtf": str(gtf),
            "reference_fai": str(fai),
            "tss_half_widths_bp": [500, 1500, 2000],
            "atac_peak_files": {"0h": "peaks/ATAC_0h_peaks.narrowPeak", "6h": "peaks/ATAC_6h_peaks.narrowPeak"},
            "histone_peak_files": {
                "H3K27ac_0h": "06.peak/h3ac0.bed",
                "H3K27ac_6h": "06.peak/h3ac6.bed",
                "H3K4me3_0h": "06.peak/h3k40.bed",
                "H3K4me3_6h": "06.peak/h3k46.bed",
            },
        }
        config_path = root / "config.json"
        write(config_path, json.dumps(config))
        environment = os.environ.copy()
        environment["SLURM_JOB_ID"] = "999"
        subprocess.run(
            [sys.executable, str(Path(__file__).with_name("export_and_count_definitions.py")), "--config", str(config_path)],
            check=True,
            env=environment,
            capture_output=True,
            text=True,
        )
        result = output_root / "export_999"
        qc = list(csv.DictReader((result / "05_v27_reproduction_QC.tsv").open(encoding="utf-8"), delimiter="\t"))
        assert {row["status"] for row in qc} == {"PASS"}, qc
        promoter_bed = result / "BED/promoter/TSS_pm2000bp/P4_four_way_common_intersection/promoter_0h.bed"
        enhancer_bed = result / "BED/enhancer/TSS_pm2000bp/H3K4me3_either_time_v27/enhancer_0h.bed"
        enhancer_consensus = result / "BED/enhancer/TSS_pm2000bp/H3K4me3_either_time_v27/enhancer_consensus_0h_6h.bed"
        assert promoter_bed.read_text(encoding="utf-8").count("\n") == 1
        assert enhancer_bed.read_text(encoding="utf-8").count("\n") == 1
        assert enhancer_consensus.read_text(encoding="utf-8").split("\t")[:3] == ["chr1", "8000", "8300"]
        assert (result / "说明_定义与输出列.txt").is_file()
        assert (result / "completed.json").is_file()
    print("INTEGRATION_TEST_OK")


if __name__ == "__main__":
    main()
