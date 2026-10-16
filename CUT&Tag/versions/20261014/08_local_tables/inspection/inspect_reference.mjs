import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const inputPath = "F:/2026.09.24_CUT&TAG/结果汇总整理/参考格式-项目数据进度总表20210715更新.xlsx";
const previewDir = "work/reference_preview";

await fs.mkdir(previewDir, { recursive: true });
const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(inputPath));

const sheetSummary = await workbook.inspect({
  kind: "sheet",
  include: "id,name",
  maxChars: 12000,
});
console.log("=== SHEETS ===");
console.log(sheetSummary.ndjson);

const broad = await workbook.inspect({
  kind: "workbook,sheet,table,drawing",
  maxChars: 20000,
  tableMaxRows: 10,
  tableMaxCols: 18,
  tableMaxCellChars: 100,
});
console.log("=== OVERVIEW ===");
console.log(broad.ndjson);

const logicalRanges = ["A1:Y94", "A1:Y67", "A1:AB48", "A1:X7"];

for (let i = 0; i < workbook.worksheets.items.length; i += 1) {
  const sheet = workbook.worksheets.getItemAt(i);
  const address = logicalRanges[i] ?? "A1:Z100";
  console.log(`=== USED ${i + 1}: ${sheet.name} ${address} ===`);
  const table = await workbook.inspect({
    kind: "table",
    sheetId: sheet.name,
    range: address,
    include: "values,formulas",
    tableMaxRows: 35,
    tableMaxCols: 22,
    tableMaxCellChars: 100,
    maxChars: 12000,
  });
  console.log(table.ndjson);

  const style = await workbook.inspect({
    kind: "computedStyle",
    sheetId: sheet.name,
    range: address,
    maxChars: 5000,
  });
  console.log(`=== STYLE ${sheet.name} ===`);
  console.log(style.ndjson);

  const preview = await workbook.render({
    sheetName: sheet.name,
    range: address,
    scale: 1,
    format: "png",
  });
  const safeName = `${String(i + 1).padStart(2, "0")}_${sheet.name.replace(/[\\/:*?"<>|]/g, "_")}.png`;
  await fs.writeFile(`${previewDir}/${safeName}`, new Uint8Array(await preview.arrayBuffer()));
}
