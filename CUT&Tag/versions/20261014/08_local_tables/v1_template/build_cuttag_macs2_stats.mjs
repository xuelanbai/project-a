import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const outputDir = "outputs/01a0d290-1cba-7113-887a-16c91462a481";
const outputPath = `${outputDir}/3D4_CUTTag_MACS2_项目数据统计表.xlsx`;
const previewPath = "work/3D4_CUTTag_MACS2_项目数据统计表.png";

const fontName = "Microsoft YaHei";
const colors = {
  navy: "#1F4E78",
  blue: "#2F75B5",
  blueLight: "#D9EAF7",
  green: "#548235",
  greenLight: "#E2F0D9",
  orange: "#C65911",
  orangeLight: "#FCE4D6",
  purple: "#7030A0",
  purpleLight: "#E4DFEC",
  gray: "#7F8C8D",
  grayLight: "#E7E6E6",
  yellow: "#FFF2CC",
  title: "#17365D",
  white: "#FFFFFF",
  text: "#1F1F1F",
  border: "#B7C9D6",
};

const workbook = Workbook.create();
const sheet = workbook.worksheets.add("项目统计");
sheet.showGridLines = false;
sheet.tabColor = colors.navy;

// 标题与项目信息
sheet.mergeCells("A1:AD1");
sheet.getRange("A1").values = [["3D4细胞 CUT&Tag 项目数据进度总表"]];
sheet.getRange("A1:AD1").format = {
  fill: colors.title,
  font: { name: fontName, size: 18, bold: true, color: colors.white },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};
sheet.getRange("A1:AD1").format.rowHeight = 32;

sheet.mergeCells("A2:AD2");
sheet.getRange("A2").values = [["主要峰识别结果：MACS2；119和137按当前样本命名保留"]];
sheet.getRange("A2:AD2").format = {
  fill: colors.blueLight,
  font: { name: fontName, size: 10, color: colors.title, italic: true },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};
sheet.getRange("A2:AD2").format.rowHeight = 22;

const infoLabels = ["A3", "G3", "K3", "O3", "S3", "W3", "A4", "G4", "K4", "O4", "S4", "W4"];
const infoLabelValues = [
  "项目", "物种", "参考基因组", "主要峰识别", "归一化", "分析流程",
  "样本类型", "状态/时间点", "测序批次", "结果状态", "统计日期", "填写说明",
];
infoLabels.forEach((addr, i) => {
  sheet.getRange(addr).values = [[infoLabelValues[i]]];
  sheet.getRange(addr).format = {
    fill: colors.navy,
    font: { name: fontName, size: 10, bold: true, color: colors.white },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  };
});

const infoMerges = [
  ["B3:F3", "3D4细胞 CUT&Tag"],
  ["H3:J3", "猪（Sus scrofa）"],
  ["L3:N3", "Sscrofa11.1"],
  ["P3:R3", "MACS2"],
  ["T3:V3", "CPM"],
  ["X3:AD3", "nf-core/cutandrun 3.2.2"],
  ["B4:F4", "3D4细胞"],
  ["H4:J4", "119、137（具体生物学含义待确认）"],
  ["L4:N4", "共3批；各样本批次待填写"],
  ["P4:R4", "MACS2峰数与TSS指标待填写"],
  ["T4:V4", new Date("2026-09-24T00:00:00+08:00")],
  ["X4:AD4", "黄色单元格为待填写；灰色表示不适用"],
];
for (const [range, value] of infoMerges) {
  sheet.mergeCells(range);
  sheet.getRange(range).values = [[value]];
  sheet.getRange(range).format = {
    fill: colors.white,
    font: { name: fontName, size: 10, color: colors.text },
    horizontalAlignment: "left",
    verticalAlignment: "center",
    wrapText: true,
  };
}
sheet.getRange("T4:V4").format.numberFormat = "yyyy-mm-dd";
sheet.getRange("A3:AD4").format.borders = { preset: "all", style: "thin", color: colors.border };
sheet.getRange("A3:AD4").format.rowHeight = 26;

sheet.mergeCells("A6:AD6");
sheet.getRange("A6").values = [["填表规则：IgG作为阴性对照，不纳入靶标峰数汇总；主统计口径为MACS2过滤后峰数。"]];
sheet.getRange("A6:AD6").format = {
  fill: colors.yellow,
  font: { name: fontName, size: 10, bold: true, color: "#7F6000" },
  horizontalAlignment: "left",
  verticalAlignment: "center",
};
sheet.getRange("A6:AD6").format.rowHeight = 24;

// 分区表头
const sectionHeaders = [
  ["A7:F7", "样本信息", colors.navy],
  ["G7:J7", "测序数据", colors.blue],
  ["K7:R7", "比对与去重复", colors.green],
  ["S7:V7", "MACS2峰统计", colors.orange],
  ["W7:Y7", "TSS富集与结论", colors.purple],
  ["Z7:AD7", "文件与备注", colors.gray],
];
for (const [range, label, fill] of sectionHeaders) {
  sheet.mergeCells(range);
  sheet.getRange(range).values = [[label]];
  sheet.getRange(range).format = {
    fill,
    font: { name: fontName, size: 10, bold: true, color: colors.white },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  };
}
sheet.getRange("A7:AD7").format.borders = { preset: "all", style: "thin", color: colors.white };
sheet.getRange("A7:AD7").format.rowHeight = 23;

const headers = [[
  "送测编号", "测序批次", "状态/时间点", "样品", "靶标/对照", "生物学重复",
  "原始数据量（Gb）", "总读段数", "R1读段数", "R2读段数",
  "比对读段数", "比对率（%）", "正确配对数", "正确配对率（%）", "重复读段数", "重复率（%）", "非重复率（%）", "有效片段数",
  "MACS2峰型", "MACS2原始峰数", "MACS2过滤后峰数", "FRiP（%）",
  "TSS中心/侧翼信号比", "TSS富集形态", "数据质量结论",
  "原始数据目录", "比对BAM路径", "MACS2峰文件路径", "TSS图路径", "备注",
]];
sheet.getRange("A8:AD8").values = headers;
const headerGroups = [
  ["A8:F8", colors.blueLight],
  ["G8:J8", "#DDEBF7"],
  ["K8:R8", colors.greenLight],
  ["S8:V8", colors.orangeLight],
  ["W8:Y8", colors.purpleLight],
  ["Z8:AD8", colors.grayLight],
];
for (const [range, fill] of headerGroups) {
  sheet.getRange(range).format = {
    fill,
    font: { name: fontName, size: 9, bold: true, color: colors.text },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
  };
}
sheet.getRange("A8:AD8").format.borders = { preset: "all", style: "thin", color: colors.border };
sheet.getRange("A8:AD8").format.rowHeight = 42;

const samples = [
  ["S01", null, "119", "3D4-119-IgG", "IgG（对照）", "—"],
  ["S02", null, "119", "3D4-119-Pol2-rep1", "Pol2", 1],
  ["S03", null, "119", "3D4-119-Pol2-rep2", "Pol2", 2],
  ["S04", null, "119", "3D4-119-Tag-rep1", "Tag", 1],
  ["S05", null, "119", "3D4-119-Tag-rep2", "Tag", 2],
  ["S06", null, "137", "3D4-137-IgG", "IgG（对照）", "—"],
  ["S07", null, "137", "3D4-137-Pol2-rep1", "Pol2", 1],
  ["S08", null, "137", "3D4-137-Pol2-rep2", "Pol2", 2],
  ["S09", null, "137", "3D4-137-Tag-rep1", "Tag", 1],
  ["S10", null, "137", "3D4-137-Tag-rep2", "Tag", 2],
];

const rows = samples.map((s) => {
  const isIgG = String(s[4]).startsWith("IgG");
  return [
    ...s,
    null, null, null, null,
    null, null, null, null, null, null, null, null,
    isIgG ? "不适用" : null,
    isIgG ? "不适用" : null,
    isIgG ? "不适用" : null,
    isIgG ? "不适用" : null,
    isIgG ? "不适用" : null,
    isIgG ? "不适用" : "待评估",
    "待评估",
    null, null, null, null, null,
  ];
});
sheet.getRange("A9:AD18").values = rows;
sheet.getRange("A9:AD18").format = {
  font: { name: fontName, size: 9, color: colors.text },
  verticalAlignment: "center",
  wrapText: false,
};
sheet.getRange("A9:AD18").format.borders = { preset: "all", style: "thin", color: colors.border };
sheet.getRange("A9:AD18").format.rowHeight = 22;
sheet.getRange("A9:C18").format.horizontalAlignment = "center";
sheet.getRange("E9:F18").format.horizontalAlignment = "center";
sheet.getRange("G9:Y18").format.horizontalAlignment = "center";

// 待填输入与计算列
for (const range of ["B9:B18", "G9:P18", "R9:R18", "S10:V13", "S15:V18", "W10:Y13", "W15:Y18", "Z9:AD18"]) {
  sheet.getRange(range).format.fill = colors.yellow;
}
sheet.getRange("Q9").formulas = [["=IF(ISNUMBER(P9),100-P9,\"\")"]];
sheet.getRange("Q9:Q18").fillDown();
sheet.getRange("Q9:Q18").format.fill = colors.greenLight;

for (const row of [9, 14]) {
  sheet.getRange(`S${row}:X${row}`).format.fill = colors.grayLight;
  sheet.getRange(`Y${row}`).format.fill = colors.yellow;
}

sheet.getRange("B9:B18").dataValidation = { rule: { type: "list", values: ["第1批", "第2批", "第3批"] } };
sheet.getRange("S10:S13").dataValidation = { rule: { type: "list", values: ["窄峰（narrowPeak）", "宽峰（broadPeak）"] } };
sheet.getRange("S15:S18").dataValidation = { rule: { type: "list", values: ["窄峰（narrowPeak）", "宽峰（broadPeak）"] } };
sheet.getRange("X10:X13").dataValidation = { rule: { type: "list", values: ["明显", "中等", "不明显", "待评估"] } };
sheet.getRange("X15:X18").dataValidation = { rule: { type: "list", values: ["明显", "中等", "不明显", "待评估"] } };
sheet.getRange("Y9:Y18").dataValidation = { rule: { type: "list", values: ["正常", "需复核", "不适用", "待评估"] } };

sheet.getRange("G9:G18").format.numberFormat = "0.0";
sheet.getRange("H9:K18").format.numberFormat = "#,##0";
sheet.getRange("L9:L18").format.numberFormat = "0.0";
sheet.getRange("M9:M18").format.numberFormat = "#,##0";
sheet.getRange("N9:N18").format.numberFormat = "0.0";
sheet.getRange("O9:O18").format.numberFormat = "#,##0";
sheet.getRange("P9:Q18").format.numberFormat = "0.0";
sheet.getRange("R9:R18").format.numberFormat = "#,##0";
sheet.getRange("T9:U18").format.numberFormat = "#,##0";
sheet.getRange("V9:V18").format.numberFormat = "0.00";
sheet.getRange("W9:W18").format.numberFormat = "0.00";

// MACS2峰数与TSS汇总
sheet.mergeCells("A21:L21");
sheet.getRange("A21").values = [["MACS2峰数与TSS富集汇总"]];
sheet.getRange("A21:L21").format = {
  fill: colors.title,
  font: { name: fontName, size: 12, bold: true, color: colors.white },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};
sheet.getRange("A21:L21").format.rowHeight = 26;

sheet.getRange("A22:L22").values = [[
  "组别", "状态/时间点", "靶标", "R1过滤后峰数", "R2过滤后峰数", "平均峰数",
  "重复间差异（%）", "共识峰数", "R1 TSS中心/侧翼比", "R2 TSS中心/侧翼比", "TSS均值", "结论",
]];
sheet.getRange("A22:L22").format = {
  fill: colors.blue,
  font: { name: fontName, size: 9, bold: true, color: colors.white },
  horizontalAlignment: "center",
  verticalAlignment: "center",
  wrapText: true,
};
sheet.getRange("A22:L22").format.rowHeight = 38;

sheet.getRange("A23:C26").values = [
  ["119-Pol2", "119", "Pol2"],
  ["119-Tag", "119", "Tag"],
  ["137-Pol2", "137", "Pol2"],
  ["137-Tag", "137", "Tag"],
];
const sourceRows = [[10, 11], [12, 13], [15, 16], [17, 18]];
for (let i = 0; i < sourceRows.length; i += 1) {
  const row = 23 + i;
  const [r1, r2] = sourceRows[i];
  sheet.getRange(`D${row}:E${row}`).formulas = [[
    `=IF(U${r1}=\"\",\"\",U${r1})`,
    `=IF(U${r2}=\"\",\"\",U${r2})`,
  ]];
  sheet.getRange(`F${row}`).formulas = [[`=IF(COUNT(D${row}:E${row})=2,AVERAGE(D${row}:E${row}),\"\")`]];
  sheet.getRange(`G${row}`).formulas = [[`=IF(COUNT(D${row}:E${row})=2,IF(AVERAGE(D${row}:E${row})=0,\"\",ABS(D${row}-E${row})/AVERAGE(D${row}:E${row})*100),\"\")`]];
  sheet.getRange(`I${row}:J${row}`).formulas = [[
    `=IF(W${r1}=\"\",\"\",W${r1})`,
    `=IF(W${r2}=\"\",\"\",W${r2})`,
  ]];
  sheet.getRange(`K${row}`).formulas = [[`=IF(COUNT(I${row}:J${row})=2,AVERAGE(I${row}:J${row}),\"\")`]];
}
sheet.getRange("L23:L26").values = [["待评估"], ["待评估"], ["待评估"], ["待评估"]];
sheet.getRange("H23:H26").format.fill = colors.yellow;
sheet.getRange("L23:L26").format.fill = colors.yellow;
sheet.getRange("L23:L26").dataValidation = { rule: { type: "list", values: ["正常", "需复核", "待评估"] } };
sheet.getRange("A23:L26").format = {
  font: { name: fontName, size: 9, color: colors.text },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};
sheet.getRange("A22:L26").format.borders = { preset: "all", style: "thin", color: colors.border };
sheet.getRange("D23:F26").format.numberFormat = "#,##0";
sheet.getRange("G23:G26").format.numberFormat = "0.0";
sheet.getRange("H23:H26").format.numberFormat = "#,##0";
sheet.getRange("I23:K26").format.numberFormat = "0.00";
sheet.getRange("A23:L26").format.rowHeight = 22;

// 方法说明
sheet.mergeCells("A29:AD29");
sheet.getRange("A29").values = [["统计口径与判读说明"]];
sheet.getRange("A29:AD29").format = {
  fill: colors.gray,
  font: { name: fontName, size: 11, bold: true, color: colors.white },
  horizontalAlignment: "left",
  verticalAlignment: "center",
};

const notes = [
  "1. 峰数以MACS2过滤后峰数为主；原始峰数、峰型（窄峰/宽峰）和FRiP同时保留，避免只看一个数字。",
  "2. MACS2与SEACR采用不同算法假设。SEACR更贴近低背景、稀疏信号的CUT&Tag/CUT&RUN；MACS2更通用、便于与ChIP-seq及既往研究比较。两种软件并非必须同时作为主结果，本表只汇总MACS2。",
  "3. TSS判读：Pol2通常应在TSS附近出现清晰中心富集；IgG应整体低且较平坦。Tag靶标含义尚未明确，需结合抗体生物学预期判断，不设置统一硬阈值。",
  "4. 119和137暂按两个状态/时间点记录。黄色单元格用于粘贴结果；共识峰数需使用同一MACS2参数、同一峰型和一致的重复合并规则。",
];
for (let i = 0; i < notes.length; i += 1) {
  const row = 30 + i;
  sheet.mergeCells(`A${row}:AD${row}`);
  sheet.getRange(`A${row}`).values = [[notes[i]]];
  sheet.getRange(`A${row}:AD${row}`).format = {
    fill: i % 2 === 0 ? colors.white : "#F6F6F6",
    font: { name: fontName, size: 9, color: colors.text },
    horizontalAlignment: "left",
    verticalAlignment: "center",
    wrapText: true,
  };
  sheet.getRange(`A${row}:AD${row}`).format.rowHeight = i === 1 || i === 2 ? 38 : 28;
}
sheet.getRange("A29:AD33").format.borders = { preset: "all", style: "thin", color: colors.border };

// 列宽与冻结窗口
const widths = {
  A: 10, B: 10, C: 12, D: 25, E: 13, F: 11,
  G: 14, H: 13, I: 12, J: 12,
  K: 13, L: 12, M: 13, N: 14, O: 13, P: 12, Q: 13, R: 13,
  S: 18, T: 15, U: 17, V: 11,
  W: 18, X: 13, Y: 14,
  Z: 24, AA: 24, AB: 26, AC: 22, AD: 24,
};
for (const [col, width] of Object.entries(widths)) {
  sheet.getRange(`${col}:${col}`).format.columnWidth = width;
}
sheet.freezePanes.freezeRows(8);
sheet.freezePanes.freezeColumns(4);

workbook.recalculate();

const check = await workbook.inspect({
  kind: "table",
  sheetId: "项目统计",
  range: "A1:AD33",
  include: "values,formulas",
  tableMaxRows: 33,
  tableMaxCols: 30,
  tableMaxCellChars: 80,
  maxChars: 20000,
});
console.log("=== KEY RANGE ===");
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
await fs.mkdir("work", { recursive: true });
const preview = await workbook.render({
  sheetName: "项目统计",
  range: "A1:AD33",
  scale: 0.8,
  format: "png",
});
await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(`OUTPUT=${outputPath}`);
console.log(`PREVIEW=${previewPath}`);
