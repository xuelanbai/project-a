import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const source = "F:/2026.09.24_CUT&TAG/结果汇总整理/参考格式-项目数据进度总表20210715更新.xlsx";
const outDir = "./reference_preview";
await fs.mkdir(outDir, { recursive: true });

const input = await FileBlob.load(source);
const workbook = await SpreadsheetFile.importXlsx(input);

const sheets = await workbook.inspect({ kind: "sheet", include: "id,name", maxChars: 5000 });
console.log("SHEETS");
console.log(sheets.ndjson);

const style = await workbook.inspect({
  kind: "computedStyle",
  sheetId: "ASFV",
  range: "A1:Y12",
  maxChars: 10000,
});
console.log("STYLE");
console.log(style.ndjson);

const preview = await workbook.render({ sheetName: "ASFV", range: "A1:Y20", scale: 1, format: "png" });
await fs.writeFile(`${outDir}/ASFV_A1_Y20.png`, new Uint8Array(await preview.arrayBuffer()));
