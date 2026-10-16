import fs from 'node:fs/promises';
import path from 'node:path';
import { gunzipSync } from 'node:zlib';
import { FileBlob, SpreadsheetFile } from '@oai/artifact-tool';

const root = 'C:/Users/Luca/Documents/Codex/2026-09-17/cut-tag';
const outputDir = path.join(root, 'outputs/01a0d290-1cba-7113-887a-16c91462a481');
const fillV2 = process.argv.includes('--fill-v2');
const source = path.join(outputDir, fillV2 ? '工作簿1_CUTTag新版参数备注_待填新版统计.xlsx' : '工作簿1_已填CUTTag统计.xlsx');
const wb = await SpreadsheetFile.importXlsx(await FileBlob.load(source));
if (fillV2) {
  await populateV2(wb, JSON.parse(gunzipSync(await fs.readFile('D:/a/cuttag_v2_table_inputs.json.gz')).toString('utf8')));
  process.exit(0);
}
if (process.argv.includes('--inspect')) {
  for (const name of ['Sheet2', 'Sheet1']) {
    const s = wb.worksheets.getItem(name);
    console.log(JSON.stringify({sheet:name, values:s.getRange(name === 'Sheet2' ? 'A1:AK14' : 'A1:L15').values}));
  }
  const png = await wb.render({sheetName:'Sheet2', range:'I1:W7', scale:1.3, format:'png'});
  await fs.writeFile(path.join(outputDir,'cuttag_source_preview.png'),new Uint8Array(await png.arrayBuffer()));
  process.exit(0);
}
const oldSheets = ['Sheet2','Sheet1'];
const before = oldSheets.map(name => {
  const s = wb.worksheets.getItem(name);
  const area = name==='Sheet2' ? 'A2:AK14' : 'A3:L15';
  return {name,area,values:JSON.stringify(s.getRange(area).values),formulas:JSON.stringify(s.getRange(area).formulas)};
});
wb.worksheets.getItem('Sheet2').getRange('A1').values = [['CUT&Tag 旧流程统计（新版数值待采集）']];
wb.worksheets.getItem('Sheet1').getRange('B1').values = [['旧流程 peak 数量对照']];
const notes = wb.worksheets.add('新版参数说明');
notes.showGridLines = false;
notes.getRange('A1').values = [['CUT&Tag 新版流程参数修改记录']];
notes.getRange('A2').values = [['主流程作业44536已成功。新版各样本统计尚待服务器导出，Sheet2和Sheet1保留旧流程数值作对照。']];
notes.getRange('A3').values = [['核实状态区分：日志已确认、提交脚本设置、待实际运行参数核对。参考基因组总长2,501,912,388 bp。']];
const usage='https://nf-co.re/cutandrun/3.2.2/docs/usage/';
const params='https://nf-co.re/cutandrun/3.2.2/parameters';
const seacr='https://github.com/FredHutch/SEACR';
const effective='https://deeptools.readthedocs.io/en/3.5.2/content/feature/effectiveGenomeSize.html';
const intersect='https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html';
const macs='https://macs3-project.github.io/MACS/docs/bdgcmp.html';
const tss='https://deeptools.readthedocs.io/en/3.5.2/content/tools/computeMatrix.html';
const rows = [
 ['参考基因组与注释','Sscrofa11.1；Ensembl 115','保持Sscrofa11.1；Ensembl 115','提交脚本设置，待运行参数核对','保证BAM、FASTA、Bowtie2索引及GTF使用同一组装与染色体名称。','猪样本使用其自身参考文件。该项为一致性核对，不是更换物种。','仍需核对实际索引路径，不能只根据FASTA文件名判断。',usage],
 ['macs_gsize（bp）',2700000000,2472047704,'旧命令及新版日志均已确认','旧值来自流程GRCh38默认值，且大于本猪参考总长。新版按实际FASTA统计非N碱基。','日志pig-v2-gsize-44536.txt：non_N_bases=2472047704；fasta_total=2501912388。','非N长度为近似值。有MAPQ过滤时，读长相关的唯一可比对长度更合适；并非精确可比对长度。',effective],
 ['normalisation_mode','CPM（原方案）','CPM（保持）','提交脚本设置，待运行参数核对','用户未确认加入外源spike-in，本项目采用按测序深度归一化。','nf-core提供CPM作为非spike-in方案。这里不是从已证实的Spikein运行改成CPM。','CPM不能校正全局靶标丰度变化；不同条件的差异还需要重复层面的计数分析。',usage],
 ['normalisation_binsize（bp）','旧实测待核；v3.2.2默认50','1','提交脚本设置，待运行参数核对','提高覆盖轨道位置分辨率，便于TSS与基因区域展示。','nf-core使用说明推荐非spike-in情况下CPM配合1 bp bin。','会增加运算与文件体积。1 bp覆盖轨道和TSS绘图10 bp bin是不同步骤。',usage],
 ['peakcaller顺序','原表以SEACR共识峰为主；实际顺序待核','macs2,seacr','提交脚本设置，待运行参数核对','按用户要求，以MACS2作为主要peak及下游报告来源，SEACR保留作补充。','nf-core将列表第一项作为主要caller，后续caller输出独立peak文件。','caller变化可能改变FRiP和峰数，不能据此直接判定样本质量变好或变坏。',params],
 ['IgG配对 / use_control','旧MACS2命令已使用IgG','保持IgG：119配119；137配137','旧命令已确认；新版待实际命令核对','控制各组非特异性背景。每组两个实验重复共用该编号的IgG样本。','samplesheet的control列指定对照；MACS2 --control参数建立背景模型。','IgG参与显著性建模，不是简单删掉所有与IgG重叠的区域。',usage],
 ['igg_scale_factor','旧实测待核；v3.2.2默认0.5','1.0','提交脚本设置，待运行参数核对','采用中性缩放因子，避免在未单独论证时把IgG信号乘以0.5。','这是本项目的分析选择，非猪样本特有的必改参数。','与SEACR norm联合记录。不能将该因子理解成MACS2 BAM的统一减法校正。',params],
 ['seacr_norm','旧实测待核；默认non','norm','提交脚本设置，待运行参数核对','在没有严格实验/对照标定的情况下，让SEACR进行对照到实验信号的归一化。','SEACR作者建议：除非两者已严格归一化（如spike-in），否则使用norm。','该策略可影响峰数，应记录在方法中，不能以峰越多作为正确依据。',seacr],
 ['replicate_threshold','旧实测待核；v3.2.2默认1','2','提交脚本设置，待运行参数核对','落实两个生物学重复都支持才保留的原则。','nf-core该参数规定共识峰至少由几个重复贡献。','nf-core共识区域边界与自定义严格共同覆盖区间可能不同，分别统计。',params],
 ['自定义重复overlap','此前补充脚本要求双向20%','至少1 bp；取消20%门槛','按用户要求修订；执行结果待采集','保留两个重复都出现的peak。仅相邻而无共同碱基不算重叠。','BEDTools默认intersect采用至少1 bp交集；另存双方支持的原始peak和共同覆盖区间。','共有区间数、rep1保留峰数、rep2保留峰数可能不同。不能重复相加作为总峰数。',intersect],
 ['min_peak_overlap（流程报告）','v3.2.2默认0.2；44536实际值待核','最新脚本0.000000001','当前脚本值，不能据此认定44536已使用','让流程重现性报告使用近似至少1 bp重叠规则。','该参数用于peak reproducibility报告。自定义2/2筛选由独立intersect实现。','若44536使用了早期脚本，流程旧报告可能仍为20%；以本次导出的实际参数和自定义统计为准。',params],
 ['MACS2 calling与表格过滤','q<=0.01；原表附加p/q筛选','calling保持q<=0.01；表格统计p>=5且q>=2（负log10值）','旧命令已确认；新版待XLS头部核对','保持单样本调用阈值，同时保留原表的附加筛选列便于比较。','narrowPeak第8列为-log10(P)，第9列为-log10(Q)。','附加筛选峰数与raw peaks分开；不把P值和负log10(P)方向写反。',params],
 ['IP/IgG重复读段处理','IP markdup保留；IgG dedup','保持原处理','旧命令确认；新版待BAM及命令核对','沿用nf-core对实验样本与IgG的默认处理，避免把真实富集重复一并去掉。','nf-core说明默认仅对照去重。是否进一步去除IP重复需要结合文库复杂度。','命名不能照抄参考文件中的nodup：本项目IP未声明去重。',usage],
 ['合并展示bigWig','原nf按重复输出CPM；44534另有旧gsize合并ppois','4条rep1-rep2 pooled vs同编号IgG的MACS2 ppois轨道','新后处理已修订；完成情况待导出核实','按用户要求集中展示119/137 × Pol2/Tag的合并信号。两个BAM联合输入后重新计算。','MACS bdgcmp ppois输出-log10(P)。p-value轨道不取两个重复p值的算术平均。','合并提高展示深度，但不替代重复支持筛选；轨道高度也不是差异结合的效应量。',macs],
 ['TSS展示设置','原表TSS指标尚待计算','全部注释基因TSS，±3 kb，绘图bin 10 bp','补充脚本设置；图及数值待采集','同时检查单重复/IgG的CPM曲线与合并IgG相对ppois曲线。','deepTools computeMatrix reference-point按同一套注释和窗口汇总信号。','不按已富集peak筛选基因后再声称全基因TSS富集。不自动等同ATAC标准TSS评分。',tss],
 ['FRiP与片段长度口径','原表FRiP待填；片段长度有既存值','按新版实际QC重新核对后填入','新版统计源待采集','FRiP需注明所用peak集与分母。片段长度取插入片段统计，不能用read length替代。','本次采集保留MultiQC、SAMtools/Picard及fragment报告来源。','旧fragment_length.tsv仅作旧来源保存，确认BAM来源一致后才能复用。',usage]
];
const headers=['参数 / 项目','旧值或旧状态','新版设置 / 数值','核实状态','修改原因','方法依据','解释与限制','依据链接'];
notes.getRange('A5:H5').values=[headers];
notes.getRange(`A6:H${5+rows.length}`).values=rows;
const body=notes.getRange(`A5:H${5+rows.length}`);
body.format.font={name:'Microsoft YaHei',size:11,color:'#222222'};
body.format.verticalAlignment='top';
body.format.wrapText=true;
notes.getRange('A1').format.font={name:'Microsoft YaHei',size:15,bold:true};
notes.getRange('A2:A3').format.font={name:'Microsoft YaHei',size:10,color:'#805500'};
notes.getRange('A5:H5').format.fill='#E2F0D9';
notes.getRange('A5:H5').format.font.bold=true;
notes.getRange('A5:H5').format.rowHeight=30;
notes.getRange(`A6:H${5+rows.length}`).format.rowHeight=65;
for (const [col,width] of [['A',30],['B',40],['C',42],['D',30],['E',54],['F',58],['G',58],['H',50]]) {
  notes.getRange(`${col}1:${col}${5+rows.length}`).format.columnWidth=width;
}
notes.getRange('B7:C7').setNumberFormat('#,##0');
notes.getRange(`D6:D${5+rows.length}`).format.fill='#FFF2CC';
notes.getRange('D7').format.fill='#FFFFFF';
notes.freezePanes.freezeRows(5);
wb.recalculate();
for (const snap of before) {
  const area=wb.worksheets.getItem(snap.name).getRange(snap.area);
  if(JSON.stringify(area.formulas)!==snap.formulas) throw new Error('Original formulas changed: '+snap.name);
  const expected=JSON.parse(snap.values), actual=area.values, formulas=area.formulas;
  for(let r=0;r<expected.length;r++) for(let c=0;c<expected[r].length;c++) {
    const a=actual[r][c],e=expected[r][c];
    if(a===e) continue;
    if(formulas[r]?.[c] && typeof a==='number' && typeof e==='number' && Math.abs(a-e)<=1e-10*Math.max(1,Math.abs(e))) continue;
    throw new Error('Original data changed: '+JSON.stringify({sheet:snap.name,r,c,expected:e,actual:a}));
  }
}
console.log((await wb.inspect({kind:'table',range:'新版参数说明!A5:D9',include:'values',tableMaxRows:5,tableMaxCols:4})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!',options:{useRegex:true,maxResults:30}})).ndjson);
for (const [name,range] of [['parameters_values','A5:D11'],['parameters_reasons','E5:H11'],['parameters_methods','A12:D21'],['parameters_limits','E12:H21']]) {
  const png=await wb.render({sheetName:'新版参数说明',range,scale:1,format:'png'});
  await fs.writeFile(path.join(outputDir,name+'.png'),new Uint8Array(await png.arrayBuffer()));
}
const out=path.join(outputDir,'工作簿1_CUTTag新版参数备注_待填新版统计.xlsx');
await (await SpreadsheetFile.exportXlsx(wb)).save(out);
console.log('EXPORTED '+out);

async function populateV2(wb, d) {
  const assert = (ok,msg) => {if(!ok) throw new Error(msg);};
  assert(d.schema==='cuttag_v2_table_inputs_v1','Unexpected input schema');
  assert(d.samples.length===10 && d.errors.length===0 && d.missing.length===0,'Incomplete source export');
  const text = k => { assert(d.text_sources[k], 'Missing source '+k); return d.text_sources[k].text; };
  const sourceKey = suffix => {const keys=Object.keys(d.text_sources).filter(k=>k.endsWith(suffix));assert(keys.length===1,'Ambiguous source '+suffix);return keys[0];};
  const tsv = s => {const lines=s.trim().split(/\r?\n/), h=lines.shift().split('\t');return lines.map(l=>Object.fromEntries(l.split('\t').map((v,i)=>[h[i],v])));};
  const index = rows => Object.fromEntries(rows.map(r=>[r.Sample,r]));
  const picard=index(tsv(text(sourceKey('/multiqc_picard_dups.txt'))));
  const fastqc=index(tsv(text(sourceKey('/multiqc_fastqc.txt'))));
  const trimmed=index(tsv(text(sourceKey('/multiqc_fastqc_1.txt'))));
  const peaks=new Map(d.peak_files.map(p=>[p.path,p]));
  const inventory=new Map(d.inventory.map(p=>[p.path,p]));
  const sampleByName=new Map(d.samples.map(s=>[s.group.replaceAll('_','-')+(s.group.endsWith('IgG')?'':'-rep'+s.replicate),s]));
  assert(sampleByName.size===10,'Duplicate sample identities');
  const sn = s => Object.fromEntries(s.split(/\r?\n/).filter(l=>l.startsWith('SN\t')).map(l=>{const f=l.split('\t');return [f[1].replace(/:$/,''),Number(f[2])];}));
  const matchNumber = (s,re) => {const m=s.match(re);assert(m,'Missing match '+re);return Number(m[1].replaceAll(',',''));};
  const repro=(group,caller)=>{const r=d.reproducibility.filter(r=>r.group===group&&r.caller===caller);assert(r.length===1,'Repro source '+group+caller);return r[0];};
  const getPeak = p => {assert(peaks.has(p),'Missing peak '+p);return peaks.get(p);};
  const s=wb.worksheets.getItem('Sheet2'), summary=wb.worksheets.getItem('Sheet1'), notes=wb.worksheets.getItem('新版参数说明');
  const preserved=['A3:D12','F3:H12','C14'].map(range=>({range,data:JSON.stringify(s.getRange(range).values)}));
  const oldRemarks=s.getRange('AJ3:AJ12').values;
  const set=(sheet,cell,value)=>{sheet.getRange(cell).clear({applyTo:'contents'});sheet.getRange(cell).formulas=[['']];sheet.getRange(cell).values=[[value]];};
  set(s,'A1','CUT&Tag 新版统计：主流程44536／合并后处理44538');
  const headers={E:'原始测序量（Gb，按FASTQ碱基数计算）',I:'平均插入片段长度（bp，markdup BAM）',J:'原始reads总数（R1+R2）',K:'FASTQ R1 reads',L:'FASTQ R2 reads',M:'比对reads（Bowtie2，过滤前）',N:'比对率（%，分母为修剪后reads）',O:'过滤后Properly paired reads（markdup）',P:'过滤后Properly paired占原始reads（%）',Q:'重复read pairs（Picard，单位：对）',R:'dup reads（%，Picard）',S:'非重复比例（1-dup rate）',V:'MACS2流程共识区间数（组，2/2）',AA:'原始GC（%，MultiQC单端汇总）',AB:'修剪移除碱基比例（%，双端含质量修剪）',AE:'MACS2共同覆盖区间数（组，≥1 bp）',AI:'已有gene热图路径（非专项TSS）',AL:'peak调用使用BAM（IgG为dedup）',AM:'本重复获另一重复支持的MACS2峰数',AN:'MACS2实际gsize（bp）',AO:'QC来源（相对results_pig_v2）',AP:'Bowtie2输入reads（修剪后R1+R2）'};
  for(const [col,h] of Object.entries(headers)) set(s,col+'2',h);
  set(s,'AQ2','Bowtie2一致比对1次（read pairs）');
  set(s,'AR2','Bowtie2一致比对多次（read pairs）');
  s.getRange('AQ2:AR12').format.font={name:'Microsoft YaHei',size:10};
  s.getRange('AQ2:AR2').format.fill='#E2F0D9';s.getRange('AQ2:AR2').format.wrapText=true;
  s.getRange('AQ2:AR12').format.columnWidth=26;s.getRange('AQ3:AR12').setNumberFormat('#,##0');
  for(const [col,h] of Object.entries({I:'平均插入片段长度\n（bp，markdup）',J:'原始总reads\n（R1+R2）',M:'Bowtie2比对reads\n（过滤前）',N:'比对率（%）\n修剪后reads为分母',O:'过滤后正确配对reads\n（markdup）',P:'过滤后配对reads\n占原始reads（%）',Q:'重复read pairs\n（Picard，单位：对）',R:'重复率（%）\nPicard',S:'非重复比例\n（1-重复率）',W:'MACS2筛选峰数\n−log10P≥5且−log10Q≥2'}))set(s,col+'2',h);
  s.getRange('AL2:AP12').format.font={name:'Microsoft YaHei',size:10};
  s.getRange('AL2:AP2').format.fill='#E2F0D9';
  s.getRange('AL2:AP2').format.wrapText=true;
  s.getRange('AL2:AL12').format.columnWidth=55;
  s.getRange('AM2:AN12').format.columnWidth=23;
  s.getRange('AO2:AO12').format.columnWidth=55;
  s.getRange('AP2:AP12').format.columnWidth=24;
  s.getRange('I2:W2').format.rowHeight=65;
  s.getRange('I2:W2').format.wrapText=true;
  s.getRange('M2:M12').format.columnWidth=24;
  const sampleChecks=[];
  for(let row=3;row<=12;row++) {
    const name=s.getRange('C'+row).values[0][0], sample=sampleByName.get(name);
    assert(sample,'Unknown workbook sample '+name);
    const id=sample.group+'_R'+sample.replicate, igg=sample.group.endsWith('_IgG');
    const statRel='02_alignment/bowtie2/target/markdup/'+id+'.stats';
    const stats=sn(text('results_pig_v2/'+statRel));
    const bamRel=statRel.replace(/\.stats$/,'.target.markdup.sorted.bam');
    assert(inventory.has(bamRel),'Missing BAM inventory '+bamRel);
    const bamCallRel=igg?bamRel.replace('/markdup/','/dedup/').replace('.markdup.','.dedup.'):bamRel;
    assert(inventory.has(bamCallRel),'Missing calling BAM '+bamCallRel);
    const trimReports=[1,2].map(n=>text('results_pig_v2/01_prealign/trimgalore/'+id+'_'+n+'.fastq.gz_trimming_report.txt'));
    const raw=trimReports.map(t=>matchNumber(t,/Total reads processed:\s+([\d,]+)/));
    const bases=trimReports.map(t=>matchNumber(t,/Total basepairs processed:\s+([\d,]+)/));
    const written=trimReports.map(t=>matchNumber(t,/Total written \(filtered\):\s+([\d,]+)/));
    assert(raw[0]===raw[1],'Unpaired raw read counts');
    const bow=text('results_pig_v2/02_alignment/bowtie2/target/log/'+id+'.bowtie2.log');
    const pairs=matchNumber(bow,/^(\d+) reads;/), aligned1=matchNumber(bow,/\s(\d+) \([^\n]+\) aligned concordantly exactly 1 time/), alignedM=matchNumber(bow,/\s(\d+) \([^\n]+\) aligned concordantly >1 times/), unaligned=matchNumber(bow,/\s(\d+) \([^\n]+\) aligned concordantly 0 times/);
    assert(pairs===aligned1+alignedM+unaligned,'Bowtie2 reconciliation '+id);
    assert(Math.abs(100*(aligned1+alignedM)/pairs-matchNumber(bow,/(\d+\.\d+)% overall alignment rate/))<.006,'Mapping rate check');
    const dup=picard[id];assert(dup,'Missing Picard '+id);
    assert(Number(dup.READ_PAIRS_EXAMINED)*2===stats['raw total sequences'],'Picard and markdup counts differ '+id);
    assert(Number(dup.READ_PAIR_DUPLICATES)*2===stats['reads duplicated'],'Duplicate count check '+id);
    set(s,'I'+row,stats['insert size average']);
    set(s,'K'+row,raw[0]);set(s,'L'+row,raw[1]);
    s.getRange('J'+row).formulas=[['=K'+row+'+L'+row]];
    set(s,'E'+row,(bases[0]+bases[1])/1e9);
    set(s,'AQ'+row,aligned1);set(s,'AR'+row,alignedM);set(s,'AP'+row,2*pairs);
    s.getRange('M'+row).formulas=[['=2*(AQ'+row+'+AR'+row+')']];
    s.getRange('N'+row).formulas=[['=M'+row+'/AP'+row+'*100']];
    set(s,'O'+row,stats['reads properly paired']);
    s.getRange('P'+row).formulas=[['=O'+row+'/J'+row+'*100']];
    set(s,'Q'+row,Number(dup.READ_PAIR_DUPLICATES));set(s,'R'+row,100*Number(dup.PERCENT_DUPLICATION));
    s.getRange('S'+row).formulas=[['=1-R'+row+'/100']];
    set(s,'AA'+row,Number(fastqc[id]['%GC']));
    set(s,'AB'+row,100*(1-(written[0]+written[1])/(bases[0]+bases[1])));
    for(const [col,n] of [['AC',1],['AD',2]]){assert(trimmed[id+'_'+n],'Missing trimmed FastQC');set(s,col+row,Number(trimmed[id+'_'+n].avg_sequence_length));}
    set(s,'X'+row,path.posix.dirname(sample.fastq_1));set(s,'Y'+row,d.results+'/'+bamRel);
    set(s,'AL'+row,d.results+'/'+bamCallRel);set(s,'AO'+row,statRel+'；01_prealign/trimgalore；04_reporting/multiqc/*_data/multiqc_picard_dups.txt');
    if(igg) {
      for(const col of ['T','U','V','W','Z','AE','AF','AG','AI','AK','AM','AN']) set(s,col+row,'—');
      set(s,'AH'+row,'未评估');
    } else {
      const macRel='03_peak_calling/04_called_peaks/macs2/'+id+'.macs2_peaks.narrowPeak';
      const mac=getPeak(macRel), sea=getPeak('03_peak_calling/04_called_peaks/seacr/'+id+'.seacr.peaks.stringent.bed');
      const con=getPeak('03_peak_calling/05_consensus_peaks/'+sample.group+'.macs2.consensus.peaks.awk.bed');
      const rep=repro(sample.group,'MACS2'), head=text('results_pig_v2/'+macRel.replace(/narrowPeak$/,'xls'));
      assert(head.includes('--gsize 2472047704')&&head.includes('--control '+sample.control+'_R1.target.dedup.sorted.bam'),'MACS2 argument mismatch '+id);
      set(s,'T'+row,sea.raw_peaks);set(s,'U'+row,mac.raw_peaks);set(s,'V'+row,con.raw_peaks);set(s,'W'+row,mac.filtered_p5_q2);
      set(s,'Z'+row,d.results+'/'+macRel);set(s,'AE'+row,rep.shared_intervals);set(s,'AM'+row,rep['rep'+sample.replicate+'_supported']);set(s,'AN'+row,2472047704);
      set(s,'AK'+row,sample.control+'_R1');
      set(s,'AF'+row,'待补：无FRiP数值');set(s,'AG'+row,'待补：无专项TSS数据');set(s,'AH'+row,'未评估，不能判定正常');
      const gene='04_reporting/deeptools_heatmaps/gene/'+id+'.plotHeatmap.pdf';
      set(s,'AI'+row,inventory.has(gene)?d.results+'/'+gene:'未收集到');
    }
    set(s,'AJ'+row,oldRemarks[row-3][0]+'；新版来源已核对。插入片段及O列统一取markdup统计；IgG调用用dedup。FRiP/TSS未提供。');
    sampleChecks.push({id,raw:raw.reduce((a,b)=>a+b),mapped:2*(aligned1+alignedM),insert:stats['insert size average']});
  }
  s.getRange('E3:E12').setNumberFormat('0.000');s.getRange('I3:I12').setNumberFormat('0.0');
  for(const col of ['J','K','L','M','O','Q','T','U','V','W','AE','AM','AN','AP'])s.getRange(col+'3:'+col+'12').setNumberFormat('#,##0');
  for(const col of ['N','P','R','AA','AB','AC','AD'])s.getRange(col+'3:'+col+'12').setNumberFormat('0.00');
  s.getRange('S3:S12').setNumberFormat('0.0000');
  set(summary,'B1','新版峰数统计：重复支持与pooled调用分列');
  set(summary,'B9','MACS2为主：共有区间不是合并BAM后的peak数');
  summary.getRange('B3:J3').clear({applyTo:'contents'});
  summary.getRange('B11:R11').clear({applyTo:'contents'});
  summary.getRange('B3:J3').values=[['组别','靶标','SEACR rep1','SEACR rep2','共同覆盖区间','rep1获支持峰','rep2获支持峰','rep1支持率','rep2支持率']];
  summary.getRange('B11:R11').values=[['组别','靶标','MACS2 rep1','MACS2 rep2','共同覆盖区间','rep1获支持峰','rep2获支持峰','rep1支持率','rep2支持率','pooled峰数','pooled p5q2峰数','nf共识区间','p5q2共同区间','合并bigWig路径','IgG对照','实际gsize','说明']];
  const groups=['3D4_119_Pol2','3D4_119_Tag','3D4_137_Pol2','3D4_137_Tag'];
  for(let i=0;i<groups.length;i++) {
    const g=groups[i],state=g.split('_')[1],mark=g.split('_')[2], mr=12+i,sr=4+i;
    for(const [caller,row] of [['SEACR',sr],['MACS2',mr]]) {
      const r=repro(g,caller);
      summary.getRange('B'+row+':H'+row).values=[[Number(state),mark,r.rep1_called,r.rep2_called,r.shared_intervals,r.rep1_supported,r.rep2_supported]];
      summary.getRange('I'+row+':J'+row).formulas=[['=G'+row+'/D'+row,'=H'+row+'/E'+row]];
      const stem='custom_tracks/reproducible_peaks_anyOverlap/'+g+'_rep1-rep2_vs_3D4_'+state+'_IgG_rep1.'+caller+'.2of2_anyOverlap.';
      for(const [suffix,key] of [['rep1_supported.bed','rep1_supported'],['rep2_supported.bed','rep2_supported'],['sharedIntervals.bed','shared_intervals']])assert(getPeak(stem+suffix).raw_peaks===r[key],'Independent peak-count disagreement');
    }
    const pooled=d.peak_files.filter(p=>p.path.startsWith('custom_tracks/macs2_ppois/'+g+'_rep1-rep2_pooled')&&p.path.endsWith('.narrowPeak'));
    assert(pooled.length===1,'Pooled peak file mismatch');
    const track=d.merged_tracks.find(t=>t.group===g);assert(track?.exists&&track.bytes>0,'Pooled BW missing');
    const phead=text('results_pig_v2/'+pooled[0].path.replace(/narrowPeak$/,'xls'));
    assert(phead.includes('-g 2472047704')&&phead.includes(g+'_R1.target.markdup.sorted.bam')&&phead.includes(g+'_R2.target.markdup.sorted.bam'),'Pooled arguments mismatch');
    summary.getRange('K'+mr+':R'+mr).values=[[pooled[0].raw_peaks,pooled[0].filtered_p5_q2,getPeak('03_peak_calling/05_consensus_peaks/'+g+'.macs2.consensus.peaks.awk.bed').raw_peaks,repro(g,'MACS2_p5q2').shared_intervals,track.path,'3D4_'+state+'_IgG_R1',2472047704,mark==='Tag'?'MACS2与SEACR差异大；须结合重复、IgG及轨道检查，不能仅凭峰数判定质量。':'pooled仅用于合并展示，重复支持以F/G/H列为准。']];
  }
  for(const range of ['B3:J7','B11:R15']){
    summary.getRange(range).format.font={name:'Microsoft YaHei',size:11};
    summary.getRange(range).format.verticalAlignment='center';
  }
  for(const range of ['B3:J3','B11:R11']){summary.getRange(range).format.fill='#E2F0D9';summary.getRange(range).format.font.bold=true;summary.getRange(range).format.wrapText=true;summary.getRange(range).format.rowHeight=40;}
  summary.getRange('B3:C15').format.columnWidth=13;
  summary.getRange('D3:N15').format.columnWidth=19;
  summary.getRange('O11:O15').format.columnWidth=55;summary.getRange('P11:Q15').format.columnWidth=25;summary.getRange('R11:R15').format.columnWidth=55;
  for(const range of ['B4:N7','B12:N15'])summary.getRange(range).format.rowHeight=27;
  for(const range of ['D4:H7','D12:H15','K12:N15','Q12:Q15'])summary.getRange(range).setNumberFormat('#,##0');
  for(const range of ['I4:J7','I12:J15'])summary.getRange(range).setNumberFormat('0.0%');
  summary.getRange('O12:R15').format.wrapText=false;
  const summaryNotes=[
    '共同覆盖区间：两个重复peak至少重叠1 bp，对交集区间去重合并。仅相邻不算重叠；一对多重叠可使区间数与获支持原始峰数不同。',
    'rep支持率=该重复获另一重复支持的原始peak数/该重复原始peak数。MACS2 p5q2指narrowPeak列8≥5且列9≥2。',
    'pooled：两个IP BAM联合输入MACS2，以同编号IgG dedup BAM为背景重新调用。pooled峰数不要求两个重复分别显著。',
    '119-Tag：rep1支持率约4.7%，rep2约2.6%；共同区间22个。137-Tag共同区间75个。暂不据此断言119特异性结合。',
    '4条pooled ppois bigWig均存在且非空，后处理44538日志完成。仅核对文件与日志，未在本地打开轨道内容。',
    '专项TSS曲线及FRiP未收集到，不能填0或判定TSS富集正常。已有gene热图不替代专项TSS中心/侧翼统计。',
    '来源：cuttag_v2_table_inputs.json.gz；采集时间2026-09-28 09:08:32（北京时间）；results_pig_v2。原流程对照仍保留在此前工作簿。'
  ];
  summary.getRange('B18:B24').values=summaryNotes.map(v=>[v]);summary.getRange('B18:B24').format.font={name:'Microsoft YaHei',size:11};summary.getRange('B18:B24').format.rowHeight=25;
  const runtime=text('results_pig_v2/pipeline_info/pipeline_report.txt');
  for(const term of ['macs_gsize: 2472047704','normalisation_binsize: 1','peakcaller: macs2,seacr','replicate_threshold: 2','min_peak_overlap: 1.0E-9','seacr_norm: norm','igg_scale_factor: 1.0'])assert(runtime.includes(term),'Unconfirmed runtime '+term);
  set(notes,'A2','新版统计已填：10个样本、4组重复支持及pooled结果。主流程44536与后处理44538完成；FRiP与专项TSS待补。');
  set(notes,'A3','核实来源：pipeline_report.txt、12份MACS2 XLS头部、QC统计、峰文件计数及44538日志。采集时间2026-09-28 09:08（北京时间）。');
  for(const r of [6,8,9,10,12,13,14,16])set(notes,'D'+r,'44536运行报告已确认');
  for(const r of [7,11,17,18])set(notes,'D'+r,'新版MACS2 XLS及运行报告已确认');
  set(notes,'D15','独立区间统计与后处理BED计数一致');
  set(notes,'D19','44538已完成，4条bigWig均存在且非空');
  set(notes,'D20','脚本已提供；包内无专项TSS结果');
  set(notes,'D21','10个markdup统计已填；FRiP缺失');
  notes.getRange('D6:D19').format.fill='#FFFFFF';
  set(notes,'C16',1e-9);notes.getRange('C16').setNumberFormat('0.00E+00');
  set(notes,'B16','v3.2.2默认0.2；旧运行实测待核');
  set(notes,'B10','旧params.yaml：seacr,macs2');
  set(notes,'G16','44536实际使用1e-9。MultiQC重现性表只含自身条目且为0，不能当作跨重复重叠率；采用独立区间重算值。');
  set(notes,'G7','MACS2单重复与pooled均确认gsize=2472047704。非N长度仍是近似有效基因组大小，并非唯一可比对长度。');
  set(notes,'C21','片段长度取新版markdup.stats；FRiP待补');
  set(notes,'G21','I列统一使用10个样本markdup BAM的SAMtools insert size average。IgG参与peak调用时另用dedup BAM，二者统计不混用。');
  notes.getRange('A24:H24').values=[['统计项','本表口径','来源','核对结果','解释','限制','处理方式','说明']];
  const defs=[
    ['原始数据量/reads','双端原始总碱基数/1e9；总reads=R1+R2','各端trimming_report.txt','10样本双端计数一致','E列为计算值，不是公司库检报告','没有库检、PCR或浓度资料','送测与实验信息保留原表空白','不根据压缩FASTQ字节数推测Gb'],
    ['Mapped reads及比对率','2×（Bowtie2 concordant一次+多次比对read pairs）','02_alignment/bowtie2/target/log/*.bowtie2.log','比对率与日志显示值相符','N列分母为修剪后输入reads（AP列）','旧表曾以原始reads×四舍五入比对率估算','本版M列改用日志整数精确计算','不能把过滤后BAM的100% mapped当作全库比对率'],
    ['片段长度与Properly paired','markdup.stats的insert size average与reads properly paired','02_alignment/bowtie2/target/markdup/*.stats','10个样本齐全','I列bp；O列reads；P列分母为原始reads','不同过滤/去重状态统计不可混用','本版统一markdup统计；调用BAM另列AL','IgG调用用dedup，不代表I列也用dedup'],
    ['重复率与重复read pairs','Q列单位对；R列Picard PERCENT_DUPLICATION×100','multiqc_picard_dups.txt及markdup.stats','重复pairs×2与重复reads一致','S列=1-R/100','R列分母为Picard examined reads，不是原始FASTQ','保留来源口径','不能把重复对数直接除以双端总reads'],
    ['修剪及GC','AB双端（原始碱基-写出碱基）/原始碱基；AA来自MultiQC raw FastQC','trimming_report及multiqc_fastqc*.txt','10样本来源齐全','AB包括接头和质量修剪，不是纯接头比例','MultiQC raw记录每样本只汇总一个端','AA明确为单端汇总；AC/AD为修剪后各端平均读长','AB取配对长度验证前的Cutadapt written计数'],
    ['MACS2共识与重复交集','V为nf共识；AE为独立交集去重区间；AM为本重复获支持原始峰','consensus.peaks.awk.bed及reproducibility','自定义8组×3种BED计数全部一致','不能互换峰边界或简单相加','一对多重叠及边界合并导致数量不同','Sheet1把pooled、nf共识、交集分开列','p5q2为附加筛选，不改变raw peaks定义'],
    ['TSS及FRiP缺失','不填0、不自动判定正常','包内inventory和text_sources','没有专项TSS输出及FRiP数值','gene热图路径已更新为新版实际存在文件','采集状态不等于TSS作业一定未提交','保留待补标记','需补充专项profile/矩阵及具有明确分母的FRiP'],
    ['样本解释','保留原表用户提供的119/137说明','原工作簿备注','未据名称重新推断生物学含义','统计阶段只按同编号IgG匹配','低峰数或大幅caller差异不能独立证明特异性','重点核查119-Tag重复一致性','后续解释需结合轨道、阴性对照和实验设计']
  ];
  notes.getRange('A25:H32').values=defs;
  notes.getRange('A24:H32').format.font={name:'Microsoft YaHei',size:11};notes.getRange('A24:H32').format.wrapText=true;notes.getRange('A24:H32').format.verticalAlignment='top';
  notes.getRange('A24:H24').format.fill='#E2F0D9';notes.getRange('A24:H24').format.font.bold=true;notes.getRange('A24:H24').format.rowHeight=30;notes.getRange('A25:H32').format.rowHeight=72;
  set(notes,'A34','统计源：D:/a/cuttag_v2_table_inputs.json.gz（用户提供）。内部源文件路径详见上述各项及Sheet2的路径列。');
  for(const p of preserved)assert(JSON.stringify(s.getRange(p.range).values)===p.data,'Unrelated original data changed: '+p.range);
  const testInput=s.getRange('AQ4').values[0][0], testMapped=s.getRange('M4').values[0][0];
  s.getRange('AQ4').values=[[testInput+1]];wb.recalculate();
  assert(s.getRange('M4').values[0][0]===testMapped+2,'Mapping formula sensitivity failed');
  s.getRange('AQ4').values=[[testInput]];
  wb.recalculate();
  for(const c of sampleChecks){const row=sampleChecks.indexOf(c)+3;assert(s.getRange('J'+row).values[0][0]===c.raw,'Raw total formula');assert(s.getRange('M'+row).values[0][0]===c.mapped,'Exact mapped count');}
  console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30}})).ndjson);
  console.log(JSON.stringify({samples:sampleChecks,groupResults:summary.getRange('B12:N15').values,unavailable:['FRiP','专项TSS']}));
  for(const [name,sheetName,range] of [['v2_samples_qc','Sheet2','I2:S12'],['v2_samples_peaks','Sheet2','T2:W12'],['v2_groups','Sheet1','B11:N15'],['v2_seacr','Sheet1','B3:J7'],['v2_verified_params','新版参数说明','A5:D21'],['v2_definitions','新版参数说明','A24:D32']]) {
    const png=await wb.render({sheetName,range,scale:1.2,format:'png'});
    await fs.writeFile(path.join(outputDir,name+'.png'),new Uint8Array(await png.arrayBuffer()));
  }
  const out=path.join(outputDir,'工作簿1_CUTTag新版统计及参数说明.xlsx');
  await (await SpreadsheetFile.exportXlsx(wb)).save(out);
  console.log('EXPORTED '+out);
}
