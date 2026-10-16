import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const inputPath = "C:/Users/Luca/Documents/Codex/2026-09-17/cut-tag/outputs/01a0d290-1cba-7113-887a-16c91462a481/工作簿1_补充片段长度与IgG对照.xlsx";
const previewPath = "work/workbook1_final_fragment_preview.png";

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(inputPath));
const sheets = await workbook.inspect({ kind: "sheet", include: "id,name", maxChars: 6000 });
console.log("=== SHEETS ===");
console.log(sheets.ndjson);

const allSheets = Array.from(workbook.worksheets);
const target = allSheets.find((s) => s.name.toLowerCase() === "sheet3") ?? allSheets[2] ?? allSheets[0];
if (!target) throw new Error("工作簿中没有工作表");

for (const range of ["H2:J12", "AN2:AP12"]) {
  const table = await workbook.inspect({ kind: "table", sheetId: target.name, range, include: "values,formulas", tableMaxRows: 11, tableMaxCols: 3, maxChars: 5000 });
  console.log(table.ndjson);
}

const preview = await workbook.render({ sheetName: target.name, range: "H1:J13", scale: 2, format: "png" });
await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));
const iggPreview = await workbook.render({ sheetName: target.name, range: "AN1:AP13", scale: 1.5, format: "png" });
await fs.writeFile("work/workbook1_final_igg_preview.png", new Uint8Array(await iggPreview.arrayBuffer()));
console.log(`TARGET_SHEET=${target.name}`);
console.log(`PREVIEW=${previewPath}`);
