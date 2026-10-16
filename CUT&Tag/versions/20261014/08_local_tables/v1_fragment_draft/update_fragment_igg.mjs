import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const inputPath = "outputs/01a0d290-1cba-7113-887a-16c91462a481/工作簿1_已填CUTTag统计.xlsx";
const fragmentPath = "fragment_length.tsv.tmp";
const outputDir = "outputs/01a0d290-1cba-7113-887a-16c91462a481";
const outputPath = `${outputDir}/工作簿1_补充片段长度与IgG对照.xlsx`;

// These two values were verified in the supplied MultiQC general-stats extract.
const insertById = new Map([
  ["3D4_119_IgG_R1", 253.6],
  ["3D4_137_IgG_R1", 260.0],
]);

const fragmentText = await fs.readFile(fragmentPath, "utf8");
const completeRows = fragmentText.trim().split(/\r?\n/).slice(1);
for (const line of completeRows) {
  const [id, raw] = line.split("\t");
  const value = Number(raw);
  if (id && Number.isFinite(value) && value > 0) insertById.set(id, value);
}

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(inputPath));
const sheet = workbook.worksheets.getItem("Sheet3");
const samples = sheet.getRange("C3:C12").values.map(([value]) => String(value ?? ""));
if (samples.length !== 10 || new Set(samples).size !== 10 || samples.some((name) => !/^3D4-(119|137)-(IgG|Pol2-rep[12]|Tag-rep[12])$/.test(name))) {
  throw new Error(`Unexpected sample list: ${samples.join(", ")}`);
}

sheet.getRange("I2").values = [["平均插入片段长度（bp）"]];
sheet.getRange("I2").format.wrapText = true;
sheet.getRange("I:I").format.columnWidth = 21;
sheet.getRange("I3:I12").format.numberFormat = "0.0";

sheet.getRange("AO2:AO15").copyFrom(sheet.getRange("AN2:AN15"), "all");
sheet.getRange("AP2:AP15").copyFrom(sheet.getRange("AN2:AN15"), "all");
sheet.getRange("AO2:AP15").clear({ applyTo: "contents" });
sheet.getRange("AO2:AP2").values = [["IgG对照样本（按组匹配）", "建议分析名称（MACS2）"]];
sheet.getRange("AO2:AP2").format.wrapText = true;
sheet.getRange("AO2:AP2").format.horizontalAlignment = "center";
sheet.getRange("AO2:AP2").format.verticalAlignment = "center";
sheet.getRange("AO2:AP2").format.font = { bold: true };
sheet.getRange("AN:AN").format.columnWidth = 57;
sheet.getRange("AO:AO").format.columnWidth = 27;
sheet.getRange("AP:AP").format.columnWidth = 48;

const lengths = [];
const labels = [];
for (const sample of samples) {
  const group = sample.match(/^3D4-(119|137)-/)[1];
  const id = sample.replaceAll("-", "_").replace(/_rep([12])$/, "_R$1") + (sample.endsWith("-IgG") ? "_R1" : "");
  lengths.push([insertById.get(id) ?? "待计算"]);
  const igg = `3D4-${group}-IgG`;
  labels.push(sample === igg
    ? ["不适用（IgG对照）", `${sample}_IgG-control`]
    : [igg, `${sample}_vs_${igg}_MACS2`]);
}
sheet.getRange("I3:I12").values = lengths;
sheet.getRange("AO3:AP12").values = labels;
sheet.getRange("AO3:AP12").format.verticalAlignment = "center";

workbook.recalculate();
for (const range of ["C2:I12", "AN2:AP12"]) {
  const view = await workbook.inspect({kind: "table", sheetId: "Sheet3", range, include: "values,formulas", tableMaxRows: 12, tableMaxCols: 7, maxChars: 12000});
  console.log(view.ndjson);
}
const errors = await workbook.inspect({kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!", options: {useRegex: true, maxResults: 100}, summary: "final formula error scan"});
console.log(errors.ndjson);

const preview = await workbook.render({sheetName: "Sheet3", range: "AN1:AP13", scale: 1.5, format: "png"});
await fs.writeFile("work/fragment_igg_preview.png", new Uint8Array(await preview.arrayBuffer()));
await fs.mkdir(outputDir, {recursive: true});
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(`OUTPUT=${outputPath}`);
