import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {Presentation, PresentationFile} from 'file:///C:/Users/Luca/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';

const workspace = 'C:/Users/Luca/Documents/Codex/2026-09-17/cut-tag';
const build = path.join(workspace, 'work/ppt_report_20261013');
const out = path.join(workspace, 'outputs/ppt_report_20261013');
const skill = 'C:/Users/Luca/.codex/plugins/cache/openai-primary-runtime/presentations/26.1007.11041/skills/presentations';
const python = 'C:/Users/Luca/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe';
process.env.RUNTIME_NODE_MODULES = 'C:/Users/Luca/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const {finalizePresentation} = await import(pathToFileURL(path.join(skill, 'container_tools/artifact_tool_utils.mjs')).href);
const font = 'Microsoft YaHei';
const pres = Presentation.create({slideSize:{width:1280,height:720}});
const tableOwners = [];
const slideCopy = [];
const source = {
  v2: path.join(workspace,'outputs/01a0d290-1cba-7113-887a-16c91462a481/工作簿1_CUTTag新版统计及参数说明.xlsx.inspect.ndjson'),
  v1: path.join(workspace,'outputs/01a0d290-1cba-7113-887a-16c91462a481/工作簿1_补充片段长度与IgG对照.xlsx'),
  q05: path.join(workspace,'work/q05_import_51897/job-51897/counts.tsv'),
  nl: path.join(workspace,'work/nolambda_import_52126/narrow_nolambda_q05-job-52126/counts.tsv'),
  qc: path.join(workspace,'outputs/hypothesis_signal_qc_v2_51892_结果解读.md'),
  screen: path.join(workspace,'work/hypothesis_screen_import/hypothesis_screen_current_peaks'),
  run: path.join(workspace,'outputs/run_cutandrun_pig_v2.slurm'),
  pooled: path.join(workspace,'outputs/postprocess_cutandtag_pig_v2.slurm'),
  sens: path.join(workspace,'outputs/cuttag_macs2_broad_nolambda_sensitivity.slurm'),
  signal: path.join(workspace,'work/hypothesis_signal_qc_import_51892/job-51892'),
};
const urls={
  nf:'https://nf-co.re/cutandrun/3.2.2/docs/usage/',
  macs:'https://pypi.org/project/MACS2/2.2.7.1/',
  seacr:'https://github.com/FredHutch/SEACR',
  bed:'https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html',
  bw:'https://github.com/deeptools/pyBigWig#compute-summary-information-on-a-range',
};
function textbox(slide,text,x,y,w,h,size=26,bold=false,color='#222222'){
  const s=slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  s.text=text;
  s.text.style={fontSize:size,typeface:font,bold,color,alignment:'left',verticalAlignment:'top',wrap:true,insets:0};
  return s;
}
function base(title,notes=''){
  const slide=pres.slides.add();slide.background.fill='#FFFFFF';
  textbox(slide,title,54,38,1170,78,44,true,'#17324D');
  textbox(slide,String(pres.slides.items.length),1192,675,40,25,17,false,'#777777');
  slide.speakerNotes.text=notes;
  slideCopy.push({title,body:[],notes});
  return slide;
}
function note(slide,text,y=585,h=70,size=25){textbox(slide,text,54,y,1170,h,size);slideCopy.at(-1).body.push(text);}
function lines(slide,entries,y=140,gap=105,size=28){entries.forEach((t,i)=>{textbox(slide,t,54,y+i*gap,1170,gap-18,size);slideCopy.at(-1).body.push(t);});}
function table(slide,values,widths,{y=142,h=360,fontSize=24}={}){
  const t=slide.tables.add({rows:values.length,columns:values[0].length,left:54,top:y,width:1170,height:h,columnWidths:widths,values});
  t.styleOptions={headerRow:false,bandedRows:false};
  t.borders.assign({style:'solid',fill:'#B5BBC2',width:0.8});
  t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({fill:'#FFFFFF',textStyle:{fontSize,typeface:font,color:'#222222'},margins:{left:10,right:10,top:6,bottom:6},anchor:'center'});
  for(let r=0;r<values.length;r++){
    t.rows[r].height=h/values.length;
    for(let c=0;c<values[0].length;c++){
      const cell=t.getCell(r,c);cell.text.style={fontSize,typeface:font,color:'#222222',alignment:c===0?'left':'center',verticalAlignment:'middle',bold:r===0};
      if(r===0)cell.fill='#E9EEF3';
    }
  }
  tableOwners.push(pres.slides.items.length);slideCopy.at(-1).body.push(values);
  return t;
}
const inspect=(await fs.readFile(source.v2,'utf8')).split(/\r?\n/).filter(Boolean).map(v=>JSON.parse(v));
const v2samples=inspect.find(v=>v.kind==='table' && v.sheet==='Sheet2').values.slice(2,12);
const v2peak=inspect.find(v=>v.kind==='table' && v.sheet==='Sheet1').values;

// 1
{
  const s=pres.slides.add();s.background.fill='#FFFFFF';
  textbox(s,'3D4 细胞 CUT&Tag 结果汇报',64,215,1152,95,58,true,'#17324D');
  textbox(s,'分析流程、参数尝试与结果解读',64,325,1152,60,34);
  textbox(s,'119 / 137 两组状态\n猪参考基因组 Sscrofa11.1\n2026-10-13',64,445,1152,140,26);
  s.speakerNotes.text='汇报截至2026-10-13已提供并核对的结果。119/137的生物学定义按用户提供的两组状态记录，具体处理与时间对应关系应以实验记录为准。';
  slideCopy.push({title:'3D4 细胞 CUT&Tag 结果汇报',body:['分析流程、参数尝试与结果解读','119 / 137 两组状态；Sscrofa11.1；2026-10-13'],notes:s.speakerNotes.text});
}
// 2
{
  const s=base('分析目的与意义','分析目的来自用户已提出的问题。互作信息作为已有研究线索记录，本汇报没有核实所提酵母双杂论文。CUT&Tag观察靶标相关的染色质信号，不能单凭峰重叠验证物理互作。');
  lines(s,[
    '定位两组 Tag 的染色质相关信号，筛选两个生物学重复共同支持的区域。',
    '结合 Pol II 与同状态 IgG，评估共占位线索及非特异性背景。',
    '检查已提出的 rRNA 相关和炎症相关关联是否得到现有数据支持。',
    '比较不同分析设置的稳定性，为结果取舍和后续验证提供依据。',
  ],145,98,28);
  note(s,'证据边界：染色质占位、基因邻近注释和蛋白物理互作属于不同证据层次。',565,85,26);
}
// 3
{
  const s=base('样本设计与对照',`来源：用户样本表与${source.v2}。共10份双端测序文库，分3批测序。每个状态有一份IgG，两个Tag重复及两个Pol II重复共用同状态IgG。外源spike-in未确认，因此当前采用CPM。`);
  table(s,[['状态','Tag','Pol II','IgG','文库数'],['119','rep1、rep2','rep1、rep2','1份','5'],['137','rep1、rep2','rep1、rep2','1份','5']],[200,290,290,230,160],{h:205});
  lines(s,['细胞与物种：3D4，猪。参考：Sscrofa11.1，Ensembl 115 注释。','共 10 份双端文库，分 3 批测序。两靶标各有两个重复，匹配同状态 IgG。','当前按 CPM 归一化。缺少已确认的 spike-in，全局绝对占位变化无法由 CPM 判定。'],385,86,27);
}
// 4
{
  const s=base('数据处理流程',`来源：${source.run}；本地nf-core/cutandrun 3.2.2源码；${urls.nf}。基础流程使用Nextflow和预安装Conda环境，计算节点离线。具体MAPQ数值未在本汇报引用的数据包中重新核实，因此不写未确认的过滤阈值。`);
  table(s,[['步骤','处理内容','主要输出'],['原始数据','FASTQ完整性、读段质量与接头检查','校验记录、FastQC / MultiQC'],['读段修剪','Trim Galore 去接头与质量修剪','修剪后双端 FASTQ'],['比对与过滤','Bowtie2 比对猪参考，过滤并标记重复','BAM、比对率、片段与重复统计'],['信号与背景','CPM覆盖轨道，同状态IgG背景','单重复 bigWig、bedGraph'],['峰调用','MACS2 与 SEACR 分别调用','narrowPeak / BED'],['重复与展示','两重复共同区间，pooled重新计算信号','共同峰集、合并展示轨道']],[180,535,455],{h:435,fontSize:24});
  note(s,'nf-core/cutandrun 3.2.2；MACS2 2.2.7.1。后续敏感性分析沿用已完成的 BAM。',600,55,24);
}
// 5
{
  const s=base('结果文件与统计口径',`来源：${source.pooled}；${urls.bed}；${urls.macs}。sharedIntervals为双方原峰相交后合并的共同碱基区间。rep1_supported、rep2_supported保留各自获另一重复支持的完整原峰。一对多关系与区间合并导致三者计数不必相等。`);
  table(s,[['结果','含义','本项目用途'],['单重复 peak','各重复独立调用的显著区域','检查重复一致性'],['共同覆盖区间','两个重复至少有 1 bp 正长度交集','主要重复支持规则'],['pooled peak','两个实验 BAM 联合输入后重新调用','补充观察合并信号'],['CPM bigWig','按测序深度归一化的覆盖信号','观察单重复及 IgG'],['MACS2 ppois bigWig','基于背景模型的 −log10(P) 分数','展示 pooled 相对 IgG 的信号']],[290,570,310],{h:370});
  note(s,'文件名 R1 / R2 在结果中表示生物学重复。共同区间数、获支持原峰数和 pooled 峰数分别统计。',560,90,26);
}
// 6
{
  const s=base('原流程至 v2 的参数调整',`来源：${source.v2}“新版参数说明”表；${source.run}。旧实测未核实项与流程默认值分开。本页仅报告已经执行的项目设置。非N长度是有效基因组大小的近似，不能等同唯一可比对长度。`);
  table(s,[['项目','原记录 / v2 设置','调整目的'],['参考与 macs_gsize','猪参考保持；2.70×10⁹ 改为 2,472,047,704','以实际猪 FASTA 非N长度替代人默认值'],['主要 caller','seacr,macs2 改为 macs2,seacr','MACS2作为主要结果，SEACR作补充'],['归一化与 bin','CPM保持；v2覆盖 bin = 1 bp','提高轨道位置分辨率'],['SEACR背景设置','v2 igg_scale_factor = 1.0；norm','明确对照缩放与归一化策略'],['重复支持','v2 threshold = 2；交集至少 1 bp','落实两个重复共同出现原则'],['MACS2阈值','narrow；q≤0.01 保持','保留严格主分析口径']],[250,500,420],{h:400,fontSize:24});
  note(s,'旧 bin、IgG缩放、SEACR norm和重复阈值的实测值未全核实。v2是项目分析选择，并非猪样本都必须如此设置。',575,82,24);
}
// 7
{
  const rows=v2samples.map(r=>[r[2].replace('3D4-',''),(r[9]/1e6).toFixed(1),Number(r[13]).toFixed(1),Number(r[17]).toFixed(1),Number(r[8]).toFixed(1)]);
  for(const [i,state] of ['119','137'].entries()){
    const s=base(`基础数据质量：${state}组`, `来源：${source.v2} Sheet2。原始总reads=R1+R2，以百万reads显示。比对率来自Bowtie2修剪后配对输入。重复率为Picard PERCENT_DUPLICATION。插入片段均值统一来自markdup BAM SAMtools统计，IgG峰调用输入另为dedup。高重复率不能单独区分PCR重复与真实强富集。`);
    table(s,[['样本','总reads（M）','比对率（%）','重复率（%）','插入片段（bp）'],...rows.slice(i*5,i*5+5)],[295,210,210,210,245],{h:330,fontSize:25});
    note(s,'Tag 比对率较低，文库重复率较高。\n插入片段列为比对后片段长度均值，不是送测文库检测长度。',535,95,27);
  }
}
// 8
{
  const s=base('原流程与 v2 的峰数对照',`来源：${source.v1}及${source.v2}。原SEACR的组级数字为旧流程共识/合并并集，不能与v2严格2/2交集比较。本页只比较单重复原始峰数。MACS2主要阈值保持q≤0.01，v2核实gsize改变，但原始峰数相同，不能把所有结果差异归于gsize。`);
  table(s,[['组别','MACS2 原流程 / v2\nrep1，rep2','SEACR 原流程\nrep1，rep2','SEACR v2\nrep1，rep2'],['119 Pol II','16,358，16,288\n两版相同','5,419，5,614','9,069，9,651'],['119 Tag','471，812\n两版相同','1,181，1,480','1,316，2,014'],['137 Pol II','11,602，10,342\n两版相同','4,367，4,323','3,205，3,963'],['137 Tag','206，253\n两版相同','1,286，984','6,111，6,116']],[220,355,295,300],{h:340});
  note(s,'MACS2原始峰数保持相同。SEACR对背景设置变化更敏感，峰数变化不能直接解释成质量改善。旧共识/并集与新版共同区间不直接比较。',530,110,26);
}
// 9
{
  const s=base('v2 的重复支持结果',`来源：${source.v2} Sheet1及${source.screen}/feature_overlap_summary.tsv。MACS2支持率=本重复获另一重复支持的原峰数/本重复原峰数。SEACR列为相同至少1 bp交集规则的共同覆盖区间。SEACR区间中位宽度119Tag1591bp、137Tag2263.5bp，MACS2分别265.5bp与291bp。不同caller的区间宽度及模型不同。`);
  table(s,[['组别','MACS2 rep1 / rep2','MACS2共同区间','rep1 / rep2支持率','SEACR共同区间'],['119 Pol II','16,358 / 16,288','13,826','82.7% / 82.2%','8,843'],['119 Tag','471 / 812','22','4.7% / 2.6%','876'],['137 Pol II','11,602 / 10,342','9,320','78.3% / 88.6%','3,178'],['137 Tag','206 / 253','75','36.4% / 29.2%','4,858']],[205,265,230,260,210],{h:325,fontSize:24});
  note(s,'119-Tag在严格MACS2下重复支持较弱。SEACR峰通常更宽、模型不同，较多峰或较高重叠率不能直接代表更可靠。',515,110,28);
}
// 10
{
  const s=base('合并展示轨道与共同峰',`来源：${source.pooled}与${source.v2} Sheet1。44538针对四组将两个目标BAM联合输入MACS2（-t bam1 bam2），各用同状态IgG，重新生成pileup与control_lambda，再用bdgcmp ppois转换为bigWig。未将两个重复的P值或ppois算术平均。共同区间由单重复峰独立求交。不能根据文件名声明目标BAM已nodup。`);
  table(s,[['组别','两重复共同区间','pooled重新调用峰','合并展示文件'],['119 Pol II','13,826','15,045','pooled vs 119 IgG'],['119 Tag','22','14','pooled vs 119 IgG'],['137 Pol II','9,320','10,431','pooled vs 137 IgG'],['137 Tag','75','140','pooled vs 137 IgG']],[250,270,280,370],{h:310});
  lines(s,['已输出 4 条 pooled MACS2 ppois bigWig，文件名包含组别、两个重复及匹配 IgG。','合并后重新估计深度、背景与显著性，pooled峰数可以增加或减少。','展示合并轨道时，仍保留单重复与 IgG 核查。区域平均ppois也不是整个区域的P值。'],485,62,25);
}
// 11
{
  const s=base('各版本尝试及分析目的',`来源：${source.run}、${source.q05}、${source.nl}、${source.sens}。q05和nolambda仅对Tag两重复运行，沿用同状态IgG、BAMPE、gsize与输入BAM策略。broad仅记录已准备方案，截至汇报未提供输出。nomodel已在原命令中出现，未形成单独的nomodel结果版本。`);
  table(s,[['版本','主要设置','目的','结果状态'],['原流程','SEACR主caller；MACS2 q≤0.01','完成基础处理与初步峰调用','已完成'],['v2 严格分析','猪gsize、背景设置与2/2支持','统一参数与统计口径','44536 / 44538完成'],['q0.05','MACS2 narrow，保留局部背景','观察较弱候选能否重复','51897完成'],['q0.05 + nolambda','MACS2 narrow，取消局部背景','检验对背景模型的敏感程度','52126完成'],['broad方案','q≤0.05；broad-cutoff≤0.10','观察较宽信号区间','尚未提供结果']],[235,430,330,175],{h:390,fontSize:23});
  note(s,'“nomodel”已在现有命令中使用。BAMPE按实际双端片段工作，不能把nomodel当作独立放松阈值的版本。',570,85,25);
}
// 12
{
  const s=base('q≤0.05：新增候选有限，119重复支持仍弱',`来源：${source.q05}及同目录BED。q0.05 rep1/rep2均独立调用，未再次用q≤0.01筛回。共同区间至少1bp。119支持原峰74/2885与73/3766，137131/358与130/406。q0.05与严格q0.01交集计数22和75，新完全不相交候选52和56。`);
  table(s,[['Tag状态','rep1峰','rep2峰','共同区间','严格q0.01\n共同区间','新增不相交区间'],['119','2,885','3,766','74','22','52'],['137','358','406','131','75','56']],[190,185,185,190,215,205],{h:210});
  lines(s,['119：单重复峰数增加明显，但获另一重复支持的比例仅约 2.6% / 1.9%。','137：共同区间由 75 增至 131，重复支持率约 36.6% / 32.0%。','新增区间属于探索性候选。较宽阈值增加灵敏度，同时需要更严格的背景解释。'],390,82,28);
}
// 13
{
  const s=base('nolambda：峰数对局部背景模型高度敏感',`来源：${source.nl}、各重复MACS2 XLS命令与README。已核对命令保留-c同状态IgG，q=0.05，narrow，local lambda disabled。对照列为先前同q0.05保留local lambda的结果。reciprocal25为各输入峰双方至少25%覆盖的二级检查。单重复支持率119为4743/15838=29.9%、4651/18928=24.6%，137为14537/20780=70.0%、14800/20605=71.8%。${urls.macs}`);
  table(s,[['Tag状态','rep1峰','rep2峰','共同区间','保留局部背景\n共同区间','双方≥25%\n共同区间'],['119','15,838','18,928','4,805','74','4,224'],['137','20,780','20,605','15,012','131','14,144']],[190,180,180,195,220,205],{h:225});
  lines(s,['同状态IgG仍参与调用，但候选区域的局部背景校正已关闭。','119有 3,550 / 4,805，137有 6,514 / 15,012 个共同区间与严格Pol II峰相交。','新增峰数量很大，多个候选位点的原始信号不高于IgG。需保留探索性解释。'],410,77,27);
}
// Detailed hypothesis interpretation excluded from this results-only report.
if(false){
  const s=base('Pol II 共占位与背景核查',`来源：${source.qc}；${source.screen}/119_Tag_vs_119_Pol2_cooccupancy.tsv；${source.signal}/cooccupancy_matched_comparison.json。MACS2 7/22，SEACR845/876。7个候选各20个预期匹配对照实际0、0、0、1、2、7、0，全部INSUFFICIENT_MATCHES。候选先由Tag峰筛出，使用同批信号观察Tag>IgG属于选择后描述，不能算独立验证。`);
  lines(s,[
    '严格 MACS2：119-Tag 的 22 个共同区间中，7 个与119-Pol II共同峰相交。',
    'SEACR：845 / 876 个119-Tag共同区间与Pol II相交，宽峰会提高重叠率。',
    '7处候选的信号值得复核，但匹配背景区域不足，这次检验无法排除背景解释。',
    'Tag与IgG的输入处理、开放染色质偏好和单份IgG均会影响特异性判断。',
  ],145,99,28);
  note(s,'现阶段可以汇报“共占位候选”。区间重叠和同批信号回看均不足以证明直接蛋白互作。',565,80,26);
}
if(false){
  const s=base('rRNA相关区域：局部线索与参考限制',`来源：${source.signal}/rRNA.exact_signal.tsv、rRNA_and_flanks_100bp_bins.exact_signal.tsv与${source.qc}。示例为BED零基半开AEMK02000489.1:43231-45101，长度1870bp。严格137Tag共同峰只相交31bp。精确均值包含有效染色体内未记录碱基为0。rRNA注释23行（含15个5S），核tRNA注释未充分表示。nolambda在此区域两状态都出现边缘相交，137的chr14新增5S相交位点CPM2.31/2.15低于IgG5.00。`);
  textbox(s,'示例：AEMK02000489.1 上一个 1,870 bp 的 rRNA 注释片段',54,137,1170,52,27);
  table(s,[['状态','Tag rep1 CPM','Tag rep2 CPM','IgG CPM'],['119','2.19','2.84','1.58'],['137','7.16','8.01','3.26']],[220,325,325,300],{y:207,h:205});
  lines(s,['137平均信号较高，但主要信号集中于片段边缘和上游，119也有信号。','严格MACS2只触及该片段边界31 bp。nolambda新增的5S位点仍可有较高IgG背景。','标准参考与注释不能完整代表rDNA阵列，重复区域的比对限制使阴性结果难以外推。'],460,65,25);
}
if(false){
  const s=base('候选基因区域：现有信号不足以支持普遍特异富集',`来源：${source.signal}/candidate_TSS_plusminus500bp.exact_signal.tsv。数值为gene-level TSS±500bp的全区域精确均值，非全基因组TSS富集评分。不同轨道分支及duplicate策略影响严格可比性，CPM为描述性展示。POLR1D编码基因上有峰不能证明与POLR1D蛋白互作。本页不声称进行了炎症通路富集检验。`);
  table(s,[['区域与状态','Tag rep1 / rep2 CPM','匹配 IgG CPM','解读'],['HIF1A，119','1.48 / 1.89','3.47','两个重复均低于IgG'],['TRAF6，119','0.98 / 1.22','0.87','差异较小'],['TRAF6，137','2.45 / 2.87','0.57','有局部候选信号'],['POLR1D，137','1.23 / 0.84','1.16','两个重复方向不一致']],[280,300,260,330],{h:320,fontSize:24});
  note(s,'当前数据没有支持“119在这组炎症候选TSS普遍特异富集”。共同候选允许生物学交叉，也需考虑技术背景。',520,112,27);
}
// 17
{
  const s=base('结果小结与未完成项目',`来源：上述核实统计与${source.qc}。本报告按用户要求只汇报处理、版本调整与观察结果，不展开互作机制或候选基因深度分析。尚未得到全基因组专项TSS曲线/评分与FRiP。broad未提供结果。`);
  lines(s,[
    '主流程已完成，四组pooled IgG相对展示轨道已生成，统计口径与文件用途已明确。',
    '119-Tag严格重复支持较弱。q0.05增加少量共同候选，nolambda使结果大幅改变。',
    'Pol II两个重复一致性较高；Tag峰数及重复支持对调用阈值和背景模型更敏感。',
    '尚未提供：FRiP、全基因组TSS曲线与评分、broad实际结果。缺失项目不填推测值。',
  ],145,99,28);
}
// 18
{
  const s=base('方法依据与参数记录',`官方依据：${urls.nf}\n${urls.macs}\n${urls.seacr}\n${urls.bed}\n${urls.bw}\n项目来源：${source.v2}\n${source.q05}\n${source.nl}\n${source.qc}\nMACS2实际命令含BAMPE、nomodel、shift=-75、extsize=150。文档说明BAMPE使用观测双端片段且不适用非零shift。运行头部记录实际片段长度（119rep1为205bp），因此不能把extsize=150解释为实际所有片段长度。该项在正式方法报告前需核实，不能据此断言已完成纠正。`);
  lines(s,[
    'nf-core/cutandrun 3.2.2：流程、对照配对与输出口径。',
    'MACS2 2.2.7.1：BAMPE、有效基因组大小、q阈值、broad及nolambda。',
    'SEACR作者文档：背景归一化与stringent峰。BEDTools：重复峰交集。',
    '运行记录：v2 44536 / 44538；信号核查51892；q0.05 51897；nolambda 52126。',
  ],145,88,26);
  note(s,'各页备注保留数据来源与统计口径。broad仅列已准备方案；nomodel并非新增版本，不补写未运行的结果。',535,105,26);
}

await fs.mkdir(build,{recursive:true});await fs.mkdir(out,{recursive:true});
const candidate=path.join(build,'candidate.pptx');
await (await PresentationFile.exportPptx(pres)).save(candidate);
console.log(`EXPORTED draft: ${candidate}; slides=${pres.slides.items.length}`);
await fs.writeFile(path.join(build,'slide_content.json'),JSON.stringify(slideCopy,null,2),'utf8');
await fs.writeFile(path.join(build,'outline.md'),slideCopy.map((x,i)=>`## ${i+1}. ${x.title}\n\n${x.body.map(v=>typeof v==='string'?v:JSON.stringify(v)).join('\n\n')}\n`).join('\n'),'utf8');
const preview=path.join(build,'previews');await fs.mkdir(preview,{recursive:true});
for(const [i,s] of pres.slides.items.entries()){
  const png=await pres.export({slide:s,format:'png',scale:1});
  await fs.writeFile(path.join(preview,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await png.arrayBuffer()));
  console.log(`RENDERED ${i+1}/${pres.slides.items.length}`);
}
const delivery=path.join(build,'deliverables');await fs.mkdir(delivery,{recursive:true});
const final=path.join(delivery,'3D4_CUTTag_结果汇报_简洁版_20261013.pptx');
const result=await finalizePresentation({
  workspaceDir:build,candidatePath:candidate,finalPath:final,
  pythonExecutable:python,
  integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...tableOwners.flatMap(n=>['--require-native-table-slide',String(n)])],
  requiredNativeTableOwnerSlides:tableOwners,
  fontPolicy:{basis:'design',families:[font]},
  verifyArtifactToolImport:true,
  receiptPath:path.join(build,'validation_results_only.json'),
});
console.log(JSON.stringify(result));
await fs.copyFile(final,path.join(out,path.basename(final)));
console.log(`DELIVERED: ${path.join(out,path.basename(final))}`);
