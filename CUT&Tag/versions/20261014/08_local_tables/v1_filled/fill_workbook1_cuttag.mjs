import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const inputWorkbook = "C:/Users/Luca/WPSDrive/1680827061/WPS企业云盘/崖州湾国家实验室/我的企业文档/工作簿1.xlsx";
const inputStats = "F:/2026.09.24_CUT&TAG/results/cuttag_table_inputs.txt";
const outputDir = "outputs/01a0d290-1cba-7113-887a-16c91462a481";
const outputPath = `${outputDir}/工作簿1_已填CUTTag统计.xlsx`;
const previewPath = "work/workbook1_filled_sheet3.png";

const project = "/mnt/nas/swine/home/jxwu/Projects/A/08.CUT&Tag/00.26.09";
const analysis = `${project}/analysis_nfcore_cutandrun_3.2.2`;
const results = `${analysis}/results`;

const text = await fs.readFile(inputStats, "utf8");

function between(startMarker, endMarker) {
  const start = text.indexOf(startMarker);
  if (start < 0) throw new Error(`Missing marker: ${startMarker}`);
  const contentStart = start + startMarker.length;
  const end = text.indexOf(endMarker, contentStart);
  if (end < 0) throw new Error(`Missing marker: ${endMarker}`);
  return text.slice(contentStart, end).trim();
}

const generalJson = between("### DATA\n", "\n## FLAGSTAT").replace(/\bNaN\b/g, "null");
const generalBlocks = JSON.parse(generalJson);
const records = {};
for (const block of generalBlocks) {
  for (const [sample, metrics] of Object.entries(block)) {
    records[sample] = { ...(records[sample] ?? {}), ...metrics };
  }
}

function finiteNumber(value) {
  if (value === null || value === undefined || value === "") return null;
  const number = Number(value);
  return Number.isFinite(number) ? number : null;
}

function numberOrLabel(value, label = "未提供") {
  return finiteNumber(value) ?? label;
}

const flagText = between("## FLAGSTAT\n", "\n## PICARD_DUPLICATION");
const flagBlocks = flagText.split(/^### /m).filter(Boolean);
for (const block of flagBlocks) {
  const [title, ...bodyLines] = block.split(/\r?\n/);
  const sample = title.replace(/\.flagstat\s*$/, "").trim();
  const body = bodyLines.join("\n");
  const rec = records[sample] ?? (records[sample] = {});
  const capture = (pattern) => {
    const match = body.match(pattern);
    return match ? Number(match[1]) : null;
  };
  rec.flagstat_total = capture(/^(\d+) \+ \d+ in total/m);
  rec.flagstat_duplicates = capture(/^(\d+) \+ \d+ duplicates/m);
  rec.flagstat_read1 = capture(/^(\d+) \+ \d+ read1/m);
  rec.flagstat_read2 = capture(/^(\d+) \+ \d+ read2/m);
  rec.flagstat_proper = capture(/^(\d+) \+ \d+ properly paired/m);
}

const picardLines = between("## PICARD_DUPLICATION\n", "\n## PEAK_COUNTS")
  .split(/\r?\n/)
  .filter((line) => line.trim());
for (const line of picardLines.slice(1)) {
  const [sample, uniquePairs, optical, nonoptical] = line.split("\t");
  const rec = records[sample] ?? (records[sample] = {});
  rec.READ_PAIRS_UNIQUE_FROM_TABLE = Number(uniquePairs);
  rec.READ_PAIRS_DUPLICATE_OPTICAL_FROM_TABLE = Number(optical);
  rec.READ_PAIRS_DUPLICATE_NONOPTICAL_FROM_TABLE = Number(nonoptical);
}

const peakLines = between("## PEAK_COUNTS\n", "\n## SEACR_CONSENSUS_COUNTS")
  .split(/\r?\n/)
  .filter((line) => line.trim());
for (const line of peakLines.slice(1)) {
  const [sample, caller, peakType, raw, filtered, path] = line.split("\t");
  const rec = records[sample] ?? (records[sample] = {});
  const prefix = caller.toLowerCase();
  rec[`${prefix}_peak_type`] = peakType;
  rec[`${prefix}_raw_peaks`] = Number(raw);
  rec[`${prefix}_filtered_peaks`] = filtered === "n.a." ? null : Number(filtered);
  rec[`${prefix}_peak_path`] = path;
}

const consensusLines = between("## SEACR_CONSENSUS_COUNTS\n", "\n## TSS_MATRIX_HEADERS")
  .split(/\r?\n/)
  .filter((line) => line.trim());
const consensus = {};
for (const line of consensusLines.slice(1)) {
  const [group, count, path] = line.split("\t");
  consensus[group] = { count: Number(count), path };
}

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(inputWorkbook));
const sheet = workbook.worksheets.getItem("Sheet3");

const extraHeaders = [
  "GC（%）",
  "去接头碱基比例（%）",
  "过滤后平均读长R1（bp）",
  "过滤后平均读长R2（bp）",
  "唯一片段对数",
  "估计文库复杂度",
  "比对错误率（%）",
  "MACS2重复共有峰数（组）",
  "FRiP（%）",
  "TSS中心/侧翼比",
  "TSS富集判断",
  "TSS图路径",
  "指纹图路径",
  "备注",
];
const extraCols = ["AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN"];

// 延伸现有表格样式，不重排原有列。
for (const col of extraCols) {
  const target = sheet.getRange(`${col}2:${col}16`);
  target.copyFrom(sheet.getRange("Z2:Z16"), "all");
  target.clear({ applyTo: "contents" });
}
sheet.getRange("AA2:AN2").values = [extraHeaders];
sheet.getRange("AA2:AN2").format.wrapText = true;
sheet.getRange("AA2:AN2").format.horizontalAlignment = "center";
sheet.getRange("AA2:AN2").format.verticalAlignment = "center";
sheet.getRange("AA2:AN2").format.rowHeight = 40;

const samples = [
  [4, "3D4-119-IgG", "3D4_119_IgG_R1", null],
  [5, "3D4-119-Pol2-rep1", "3D4_119_Pol2_R1", "3D4_119_Pol2"],
  [6, "3D4-119-Pol2-rep2", "3D4_119_Pol2_R2", "3D4_119_Pol2"],
  [7, "3D4-119-Tag-rep1", "3D4_119_Tag_R1", "3D4_119_Tag"],
  [8, "3D4-119-Tag-rep2", "3D4_119_Tag_R2", "3D4_119_Tag"],
  [9, "3D4-137-IgG", "3D4_137_IgG_R1", null],
  [10, "3D4-137-Pol2-rep1", "3D4_137_Pol2_R1", "3D4_137_Pol2"],
  [11, "3D4-137-Pol2-rep2", "3D4_137_Pol2_R2", "3D4_137_Pol2"],
  [12, "3D4-137-Tag-rep1", "3D4_137_Tag_R1", "3D4_137_Tag"],
  [13, "3D4-137-Tag-rep2", "3D4_137_Tag_R2", "3D4_137_Tag"],
];

for (const [row, displayName, sampleId, groupId] of samples) {
  const rec = records[sampleId];
  if (!rec) throw new Error(`Missing metrics for ${sampleId}`);
  const trimmedR1 = records[`${sampleId}_1`] ?? {};
  const trimmedR2 = records[`${sampleId}_2`] ?? {};
  const isIgG = displayName.includes("IgG");
  const rawReadsPerMate = finiteNumber(rec.total_sequences);
  if (rawReadsPerMate === null) throw new Error(`Missing raw read count for ${sampleId}`);
  const totalReads = rawReadsPerMate * 2;
  const rawGb = totalReads * (finiteNumber(rec.avg_sequence_length) ?? 150) / 1e9;
  const uniquePairs = finiteNumber(rec.READS_IN_UNIQUE_PAIRS ?? rec.READS_IN_UNIQUE_PAIRS_FROM_TABLE);
  const duplicatePairs = finiteNumber(rec.READ_PAIR_DUPLICATES) ?? 0;
  const duplication = finiteNumber(rec.PERCENT_DUPLICATION);
  const tssPdf = `${results}/04_reporting/deeptools_heatmaps/gene/${sampleId}.plotHeatmap.pdf`;
  const fingerprintPdf = `${results}/04_reporting/deeptools_qc/${sampleId}.plotFingerprint.pdf`;
  const bamPath = `${results}/02_alignment/bowtie2/target/markdup/${sampleId}.target.markdup.sorted.bam`;
  const rawDir = `${project}/${displayName}`;
  const noteBase = displayName.includes("119")
    ? "119组（Pol2互作）"
    : "137组（核糖体RNA聚合酶互作，可作对照）";

  sheet.getRange(`E${row}`).values = [[rawGb]];
  sheet.getRange(`I${row}`).values = [[numberOrLabel(rec.insert_size_average)]];
  sheet.getRange(`K${row}:L${row}`).values = [[rawReadsPerMate, rawReadsPerMate]];
  sheet.getRange(`J${row}`).formulas = [[`=K${row}+L${row}`]];
  sheet.getRange(`N${row}`).values = [[numberOrLabel(rec.overall_alignment_rate)]];
  sheet.getRange(`M${row}`).formulas = [[`=ROUND(J${row}*N${row}/100,0)`]];
  sheet.getRange(`O${row}`).values = [[numberOrLabel(rec.flagstat_total)]];
  sheet.getRange(`P${row}`).formulas = [[`=IF(J${row}>0,O${row}/J${row}*100,\"\")`]];
  sheet.getRange(`Q${row}`).values = [[duplicatePairs]];
  sheet.getRange(`R${row}`).values = [[duplication === null ? "未提供" : duplication * 100]];
  sheet.getRange(`S${row}`).formulas = [[`=IF(ISNUMBER(R${row}),1-R${row}/100,\"\")`]];

  sheet.getRange(`T${row}:W${row}`).values = [[
    isIgG ? "n.a." : numberOrLabel(rec.seacr_raw_peaks),
    isIgG ? "n.a." : numberOrLabel(rec.macs2_raw_peaks),
    isIgG ? "n.a." : numberOrLabel(consensus[groupId]?.count),
    isIgG ? "n.a." : numberOrLabel(rec.macs2_filtered_peaks),
  ]];

  sheet.getRange(`X${row}:Z${row}`).values = [[
    rawDir,
    bamPath,
    isIgG ? "n.a." : rec.macs2_peak_path,
  ]];

  sheet.getRange(`AA${row}:AN${row}`).values = [[
    numberOrLabel(rec.percent_gc),
    numberOrLabel(rec.percent_trimmed),
    numberOrLabel(trimmedR1.avg_sequence_length),
    numberOrLabel(trimmedR2.avg_sequence_length),
    uniquePairs ?? "未提供",
    numberOrLabel(rec.ESTIMATED_LIBRARY_SIZE),
    finiteNumber(rec.error_rate) === null ? "未提供" : finiteNumber(rec.error_rate) * 100,
    isIgG ? "n.a." : "待计算",
    isIgG ? "n.a." : "待计算",
    isIgG ? "n.a." : "待计算",
    isIgG ? "不适用" : "待计算",
    isIgG ? "n.a." : tssPdf,
    fingerprintPdf,
    `${noteBase}${isIgG ? "；IgG阴性对照" : ""}`,
  ]];
}

// 数字格式与待补算标记。
sheet.getRange("E4:E13").format.numberFormat = "0.00\" G\"";
sheet.getRange("I4:I13").format.numberFormat = "0.0";
sheet.getRange("J4:M13").format.numberFormat = "#,##0";
sheet.getRange("N4:N13").format.numberFormat = "0.00";
sheet.getRange("O4:O13").format.numberFormat = "#,##0";
sheet.getRange("P4:P13").format.numberFormat = "0.00";
sheet.getRange("Q4:Q13").format.numberFormat = "#,##0";
sheet.getRange("R4:R13").format.numberFormat = "0.00";
sheet.getRange("S4:S13").format.numberFormat = "0.0000";
sheet.getRange("T4:W13").format.numberFormat = "#,##0";
sheet.getRange("AA4:AB13").format.numberFormat = "0.00";
sheet.getRange("AC4:AD13").format.numberFormat = "0.0";
sheet.getRange("AE4:AF13").format.numberFormat = "#,##0";
sheet.getRange("AG4:AG13").format.numberFormat = "0.000";
sheet.getRange("AA4:AK13").format.horizontalAlignment = "center";
sheet.getRange("AA4:AN13").format.verticalAlignment = "center";
sheet.getRange("AH4:AK13").format.fill = "#FFF2CC";
for (const row of [4, 9]) {
  sheet.getRange(`AH${row}:AK${row}`).format.fill = "#E7E6E6";
}

const widths = {
  AA: 10, AB: 17, AC: 19, AD: 19, AE: 15, AF: 15, AG: 14,
  AH: 20, AI: 12, AJ: 16, AK: 14, AL: 28, AM: 28, AN: 28,
};
for (const [col, width] of Object.entries(widths)) {
  sheet.getRange(`${col}:${col}`).format.columnWidth = width;
}

// 在空白行记录项目级报告路径。
sheet.getRange("A15").values = [["样本相关性图"]];
sheet.getRange("B15").values = [[`${results}/04_reporting/deeptools_qc/all_target_bams.plotCorrelation.pdf`]];
sheet.getRange("D15").values = [["PCA图"]];
sheet.getRange("E15").values = [[`${results}/04_reporting/deeptools_qc/all_target_bams.plotPCA.pdf`]];
sheet.getRange("A15").format.font = { bold: true };
sheet.getRange("D15").format.font = { bold: true };

workbook.recalculate();

const check = await workbook.inspect({
  kind: "table",
  sheetId: "Sheet3",
  range: "A1:AN16",
  include: "values,formulas",
  tableMaxRows: 16,
  tableMaxCols: 40,
  tableMaxCellChars: 90,
  maxChars: 25000,
});
console.log("=== FILLED SHEET3 ===");
console.log(check.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 300 },
  summary: "final formula error scan",
});
console.log("=== FORMULA ERRORS ===");
console.log(errors.ndjson);

await fs.mkdir(outputDir, { recursive: true });
const preview = await workbook.render({ sheetName: "Sheet3", range: "A1:AN16", scale: 0.75, format: "png" });
await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(`OUTPUT=${outputPath}`);
console.log(`PREVIEW=${previewPath}`);
