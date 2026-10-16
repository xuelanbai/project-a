import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {FileBlob,PresentationFile} from 'file:///C:/Users/Luca/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';
const root='C:/Users/Luca/Documents/Codex/2026-09-17/cut-tag';
const build=path.join(root,'work/ppt_focus_20261013');
const output=path.join(root,'outputs/ppt_report_20261013');
const skill='C:/Users/Luca/.codex/plugins/cache/openai-primary-runtime/presentations/26.1007.11041/skills/presentations';
const python='C:/Users/Luca/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe';
process.env.RUNTIME_NODE_MODULES='C:/Users/Luca/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const {finalizePresentation,applyPresentationChartFont}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const p=await PresentationFile.importPptx(await FileBlob.load('D:/Lucas临时文件/会议汇报/工作总结-10.pptx'));
const data=JSON.parse(await fs.readFile(path.join(build,'data.json'),'utf8'));
const font='Microsoft YaHei';
const blue='#4472C4',red='#C00000',gray='#7A7A7A';
const sourceBase=p.slides.items[0];
const inspected=await p.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.text','elements.position']});
const initial={records:JSON.parse(JSON.stringify(inspected.records))};
for(const r of initial.records.slice(1))p.delete(r.id);
for(const e of initial.records[0].elements)if(e.position.top>105)p.delete(e.id);
const sourceTitle=p.resolve(initial.records[0].elements.find(s=>s.text&&s.position.top<100).id);
if(!sourceTitle)throw new Error('Template heading not found');
sourceTitle.text='CUT&Tag 结果';
sourceTitle.text.style={fontSize:32,typeface:font,bold:true,color:blue,verticalAlignment:'middle',insets:0};
sourceTitle.position={left:173,top:36,width:1020,height:51};
const tables=[],charts=[],copy=[];
function text(s,t,x=84,y=125,w=1120,h=60,size=26,color='#111111',bold=false){
 const o=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 o.text=t;o.text.style={fontSize:size,typeface:font,color,bold,wrap:true,alignment:'left',verticalAlignment:'top',insets:0};return o;
}
function slide(title,notes){
 const s=sourceBase.duplicate();s.moveTo(p.slides.items.length-1);
 const titleShape=s.shapes.items.find(x=>Math.abs(x.position.left-173)<0.01&&Math.abs(x.position.top-36)<0.01);
 titleShape.text=title;
 s.speakerNotes.text=notes;
 copy.push({title,notes});
 text(s,String(copy.length),1190,685,40,22,18,'#777777');
 return s;
}
function table(s,values,widths,{x=85,y=160,h=300,fontSize=24}={}){
 const w=widths.reduce((a,b)=>a+b,0);
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,columnWidths:widths,values});
 t.styleOptions={headerRow:false,bandedRows:false};t.borders.assign({style:'solid',fill:'#CAD0D7',width:0.6});
 t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({margins:{left:9,right:9,top:4,bottom:4},anchor:'center'});
 for(let r=0;r<values.length;r++){t.rows[r].height=h/values.length;for(let c=0;c<values[0].length;c++){const cell=t.getCell(r,c);cell.fill=r===0?'#4F81BD':r%2===0?'#EAF1F6':'#FFFFFF';cell.text.style={fontSize,typeface:font,color:r===0?'#FFFFFF':'#111111',bold:r===0,alignment:c===0?'left':'center',verticalAlignment:'middle'};}}
 tables.push(copy.length);return t;
}
function chart(s,type,categories,series,{x=85,y=155,w=1100,h=370,title='',percent=false,max=undefined,stacked=false,labels=true}={}){
 const c=s.charts.add(type,{position:{left:x,top:y,width:w,height:h},title:title||undefined,titleTextStyle:{fontSize:23,bold:true},categories,series,
  hasLegend:true,legend:{position:'bottom',overlay:false,textStyle:{fontSize:21}},
  barOptions:{direction:'column',grouping:stacked?'stacked':'clustered',gapWidth:55},lineOptions:{smooth:false},
  xAxis:{textStyle:{fontSize:21},line:{fill:'#777777',width:1}},
  yAxis:{min:0,max,numberFormatCode:percent?'0%':'0.##',textStyle:{fontSize:20},majorGridlines:{fill:'#DDDDDD',width:0.5}},
  dataLabels:{showValue:labels,position:stacked?'center':'outEnd',numberFormatCode:percent?'0%':'0.##',textStyle:{fontSize:21}},
  chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF'});
 applyPresentationChartFont(c,{fontFamily:font});charts.push(copy.length);return c;
}
async function image(s,file,x,y,w,h){const b=await fs.readFile(file);s.images.add({blob:b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),contentType:'image/png',alt:path.basename(file),fit:'contain',position:{left:x,top:y,width:w,height:h}});}
const screen=path.join(root,'work/hypothesis_screen_import/hypothesis_screen_current_peaks');
const signal=path.join(root,'work/hypothesis_signal_qc_import_51892/job-51892');
const notesCommon=`数据：${screen}；${signal}。严格峰为MACS2 q≤0.01且两个重复正长度相交后的共同覆盖区间，统计单位不是合并BAM重新调用峰。注释Sscrofa11.1/Ensembl115。`;
// 1
{
 const s=slide('CUT&Tag 结果汇报：两种 Tag 的区域与基因关联',notesCommon+'研究方向来自用户先前提供的PPI候选，为预先提出的关联问题，不是本数据已验证的偏好。POLR1D背景依据https://www.ncbi.nlm.nih.gov/gene/51082/。');
 text(s,'汇报重点',85,125,1110,42,30,blue,true);
 text(s,'① 两种 Tag 的共同峰主要分布在哪里，对应哪些基因？\n② 137-Tag 与 Pol II、rRNA 注释区域是否重叠？\n③ 119-Tag 是否在炎症相关基因区域出现峰或较高信号？',85,185,1110,205,29);
 table(s,[['数据设计','当前结果口径'],['3D4猪细胞，119 / 137两组','Tag、Pol II各2个重复，同状态各1份IgG'],['主结果：MACS2严格共同峰','SEACR及放松参数结果作为补充'],['已有功能线索','来自PPI候选和基因注释，尚无正式GO/通路富集结果']],[450,660],{y:407,h:205,fontSize:23});
 text(s,'POLR1D属于Pol I/III。137–Pol II共占位与rRNA关联分别汇报。',85,632,1110,36,23,red);
}
// 2
{
 const s=slide('分析流程与版本设置',notesCommon+`实际流程与参数：${root}/outputs/run_cutandrun_pig_v2.slurm；${root}/outputs/postprocess_cutandtag_pig_v2.slurm。`);
 text(s,'FASTQ质控与修剪 → 猪参考比对及过滤 → 同状态IgG背景峰调用 → 重复交集与区域注释',85,128,1110,80,27,blue,true);
 table(s,[['版本','主要调整','目的'],['原流程','MACS2 q≤0.01，gsize=2.70×10⁹','完成基础处理'],['猪参数v2','gsize=2,472,047,704，CPM，2/2重复支持','统一猪参考和统计口径'],['q0.05','仅放宽MACS2显著性阈值','观察较弱候选'],['q0.05 + nolambda','保留IgG，取消局部lambda背景估计','查看对背景模型的敏感程度'],['broad','已准备方案，尚未提供运行结果','本汇报不列结果数值']],[235,560,315],{y:228,h:310,fontSize:23});
 text(s,'展示轨道：两个实验BAM联合输入MACS2，匹配同状态IgG，生成pooled ppois bigWig。\n区域分析：仍以单重复峰的共同区间为准。',85,573,1110,87,25);
}
// 3
{
 const s=slide('各版本共同峰数',notesCommon+'q0.05：job51897/counts.tsv；nolambda：job52126/counts.tsv。SEACR为stringent共同峰。');
 table(s,[['Tag组别','MACS2 q0.01','MACS2 q0.05','q0.05 + nolambda','SEACR'],['119','22','74','4,805','876'],['137','75','131','15,012','4,858']],[210,220,220,270,190],{y:165,h:220,fontSize:25});
 text(s,'主要观察',85,438,1110,40,28,blue,true);
 text(s,'119在严格MACS2下共同峰较少。q0.05新增52个与严格集合不相交的共同区间。\n137由75个增至131个，新增56个不相交区间。\nnolambda峰数大幅增加，说明结果对局部背景设置敏感。',85,495,1110,143,26);
}
// 4
{
 const names=['MACS2 119','MACS2 137','SEACR 119','SEACR 137'];const keys=['macs2_119_Tag','macs2_137_Tag','seacr_119_Tag','seacr_137_Tag'];
 const s=slide('共同峰主要分布在 TSS 邻近区域',notesCommon+'本页为当前BED与Ensembl115 gene-level TSS新汇总的描述性分布。优先分类：与任一TSS±2kb正长度相交、其余与基因体相交、其余区域。不可等同完整外显子/内含子注释。分布比例没有与可比对背景作富集检验。');
 const series=[['TSS ±2 kb','promoter',blue],['其余基因体','body','#86BBD8'],['其他区域','distal','#B9BEC5']].map(([name,key,fill])=>({name,fill,values:keys.map(k=>data.all[k].regions[key]/data.all[k].peaks.length)}));
 chart(s,'bar',names,series,{y:145,h:360,percent:true,max:1,stacked:true});
 table(s,[['组别','TSS±2 kb','其余基因体','其他区域'],...keys.map((k,i)=>[names[i],...Object.values(data.all[k].regions).map(String)])],[285,275,275,275],{y:519,h:143,fontSize:21.4});
}
// 5
{
 const s=slide('严格共同峰的对应基因与区域示例',notesCommon+'最近基因来自已有nearest_gene.tsv，去重按Ensembl gene ID且剔除.-1无注释值。不同峰可指向同一基因、一个峰也可能对应多个等距基因。最近基因不能直接认定调控靶基因。');
 text(s,'119：22个共同区间，对应17个不同的最近基因ID。137：75个共同区间，对应72个ID。',85,128,1110,62,25,blue,true);
 const chosen=[['119','DNAJA3'],['119','MEIOB'],['119','HSPA6'],['137','RPS29'],['137','PCSK5'],['137','VPS13A']];
 const rows=chosen.map(([st,g])=>{const r=data.all[`macs2_${st}_Tag`].near.find(a=>a[6]===g);return [st,g,`${r[0]}:${r[1]}–${r[2]}`,+r[9]===0?'与基因体相交':`${r[9]} bp`];});
 table(s,[['Tag','最近基因','共同区间（BED坐标）','距基因边界'],...rows],[115,220,485,290],{y:220,h:330,fontSize:23});
 const unplaced=k=>Object.entries(data.all[k].chr).filter(([c])=>c.startsWith('AEMK')).reduce((a,[,n])=>a+n,0);
 text(s,`未定位contig上：119有${unplaced('macs2_119_Tag')} / 22个，137有${unplaced('macs2_137_Tag')} / 75个。\n这些区域的注释和重复序列背景会影响对应基因解释。`,85,585,1110,78,25);
}
// 6
{
 const s=slide('137 与 Pol II：严格共同峰中 53 / 75 个相交',notesCommon+'各Tag共同区间对同状态PolII共同区间至少1bp正长度重叠，每个Tag区间仅计一次。PolII共占位不能验证蛋白互作。');
 const keys=['macs2_119_Pol2','macs2_137_Pol2','seacr_119_Pol2','seacr_137_Pol2'];
 chart(s,'bar',['MACS2 119','MACS2 137','SEACR 119','SEACR 137'],[{name:'Tag共同区间与同状态Pol II相交比例',values:keys.map(k=>data.compare[k].hit/data.compare[k].total),fill:blue}],{y:146,h:350,max:1,percent:true});
 table(s,[['MACS2 119','MACS2 137','SEACR 119','SEACR 137'],keys.map(k=>`${data.compare[k].hit} / ${data.compare[k].total}`)],[277.5,277.5,277.5,277.5],{y:512,h:84,fontSize:24});
 text(s,'两个Tag也有交叉：119共同峰中，MACS2有9/22、SEACR有848/876与137相交。\n这里报告区域重叠比例，尚未获得相对背景的显著性结论。',85,614,1110,70,23);
}
// 7
{
 const s=slide('137 与 rRNA：有少量注释相交，主要为局部信号',notesCommon+'rRNA相交来自feature_overlap_summary.tsv及相交BED。AEMK02000489.1坐标为零基半开。rRNA注释23条，其中15条5S，14条bigWig均不含MT。');
 table(s,[['共同峰集','峰总数','rRNA相交','5S相交','tRNA相交'],['MACS2 119-Tag','22','1','0','0'],['MACS2 137-Tag','75','1','0','0'],['SEACR 119-Tag','876','0','0','0'],['SEACR 137-Tag','4,858','2','0','0']],[310,200,200,200,200],{y:150,h:235,fontSize:24});
 text(s,'MACS2 137示例',85,419,1110,39,28,blue,true);
 text(s,'共同峰：AEMK02000489.1:42318–43262\nrRNA注释：AEMK02000489.1:43231–45101\n相交31 bp，覆盖该1,870 bp注释片段的1.7%。119也在该片段边缘相交。',85,473,1110,130,25);
 text(s,'当前注释未完整表示rDNA串联阵列，不能据此推断全部rDNA的占位情况。',85,637,1110,38,23,red);
}
// 8
{
 const s=slide('rRNA 注释片段附近的 137-Tag 信号',notesCommon+'图值直接来自rRNA_and_flanks_100bp_bins.exact_signal.tsv，100bp分箱mean_full_interval。区域AEMK02000489.1:43231-45101，图仅绘137两个重复与IgG。CPM轨道为描述性，未使用误差线或显著性标记。');
 const rows=data.rrnaBins;const cat=rows.map((r,i)=>i%5===0?String(r.offset_bp):'');
 chart(s,'line',cat,[['137 Tag rep1','137_Tag_rep1',red],['137 Tag rep2','137_Tag_rep2','#E88869'],['137 IgG','137_IgG',gray]].map(([name,key,col])=>({name,values:rows.map(r=>+r[`${key}.mean_full_interval`]),line:{fill:col,width:2.5},marker:{symbol:'none'}})),{y:160,h:345,labels:false,title:'100 bp 分箱平均 CPM（横轴为距rRNA起始处的bp数）'});
 table(s,[['整个1,870 bp区域','119 Tag rep1 / rep2','119 IgG','137 Tag rep1 / rep2','137 IgG'],['平均CPM','2.19 / 2.84','1.58','7.16 / 8.01','3.26']],[245,260,130,300,175],{y:524,h:88,fontSize:21.4});
 text(s,'较强信号集中于起始边缘及上游，区域中段接近零。137有局部较高信号，尚不足以概括完整rDNA富集。',85,636,1110,48,23);
}
// 9
{
 const s=slide('POLR1D 基因区域的 IGV 展示',notesCommon+'图片为用户提供11_5044329-5089766_strand+_ENSSSCG00000009310_POLR1D_137.png，原图完整嵌入，坐标显示chr11:5042187-5049984。图含pooled ppois与单重复CPM等不同类型轨道，不比较不同轨道类型绝对高度。文件标识POLR1D位置，Refseq下方标签为空。');
 await image(s,'D:/Lucas临时文件/会议汇报/11_5044329-5089766_strand+_ENSSSCG00000009310_POLR1D_137.png',70,128,880,541);
 text(s,'图中观察',970,145,245,45,26,blue,true);
 text(s,'Pol II轨道有局部峰。\n\nTag有较弱信号，所显示的Tag narrowPeak轨道未见相应峰。\n\n这是POLR1D编码基因区域，和rRNA区域是两个位置问题。',970,202,240,340,23);
 text(s,'图片：已有pooled与单重复轨道快照',970,589,240,65,21.4,gray);
}
// 10
{
 const s=slide('炎症相关候选：哪些基因区域出现共同峰',notesCommon+'候选基因以gene body±2kb为窗口，命中来自candidate_gene_window_hits.tsv。原先共同轴与119/137偏好分支是先验PPI线索，不是CUTTag证实的特异性。');
 table(s,[['候选窗口','MACS2 119','MACS2 137','SEACR 119','SEACR 137'],...['HIF1A','JAK2','MYD88','TRAF6','CHUK','POLR2H','POLR1D'].map(g=>[g,...['macs2_119','macs2_137','seacr_119','seacr_137'].map(k=>data.candidates[k].includes(g)?'有':'无')])],[250,215,215,215,215],{y:150,h:390,fontSize:24});
 text(s,'119的SEACR共同峰命中HIF1A、JAK2、MYD88，137也命中这些区域。\n这些窗口在严格MACS2共同峰中均未命中。当前结果允许两组共享候选。',85,576,1110,81,25);
}
// 11
{
 const s=slide('炎症候选 TSS 信号与 IgG 对照',notesCommon+'图值为candidate_TSS_plusminus500bp.exact_signal.tsv的mean_full_interval。gene-level TSS±500bp区域均值。两个状态分别用自身IgG，不把测序深度归一化等同跨状态绝对占位校准。');
 const genes=['HIF1A','JAK2','MYD88','TRAF6'];const rows=genes.map(g=>data.tss.find(r=>r.name===g));
 for(const [i,state] of ['119','137'].entries())chart(s,'bar',genes,[['Tag rep1',`${state}_Tag_rep1`,blue],['Tag rep2',`${state}_Tag_rep2`,'#8FAADC'],['IgG',`${state}_IgG`,gray]].map(([name,k,fill])=>({name,fill,values:rows.map(r=>+r[`${k}.mean_full_interval`])})),{x:75+i*605,y:150,w:550,h:395,title:`${state}：TSS ±500 bp 平均CPM`,labels:false});
 text(s,'119：HIF1A两个重复均低于IgG。TRAF6两个重复略高于IgG（0.98、1.22 vs 0.87）。\n137：JAK2、MYD88和TRAF6存在两个重复均高于IgG的局部描述性信号。',85,587,1110,86,24);
}
// 12
{
 const s=slide('TRAF6 区域的 IGV 展示',notesCommon+'用户本地TRAF6-119.png原图完整嵌入。图片显示chr2:24601253-24609585。pooled Tag/PolII和单重复PolII轨道混排，纵轴不完全一致。图中的137Tag峰在TRAF6上游，而不是同一个TSS±500bp窗口。');
 await image(s,'D:/Lucas临时文件/会议汇报/TRAF6-119.png',70,127,880,541);
 text(s,'图中观察',970,145,240,42,26,blue,true);
 text(s,'137-Tag在上游可见pooled峰，并与Pol II局部信号相邻或重叠。\n\n119-Tag有弱信号，图中未显示对应pooled peak。\n\n截图区域与上一页TSS±500 bp窗口不同。',970,202,240,355,23);
 text(s,'仅作已有区域实例，不代表全基因组功能富集。',970,588,240,72,21.4,gray);
}
// 13
{
 const s=slide('对应基因的功能关联与富集结果状态',notesCommon+'SEACR最近基因见nearest_genes/*.nearest_gene.tsv。功能方向为用户提供的候选功能分类，属于描述性归类。未执行新的GO/KEGG/Reactome富集。GO方法背景：https://www.geneontology.org/docs/go-enrichment-analysis/。现有gene_list_for_enrichment.tsv是输入列表不是富集结果。');
 const selected=['IL6','NFKB2','STAT3','TRAF2','TRAF7','HIF1A','JAK2'];
 table(s,[['功能关联线索','119 SEACR最近基因','137 SEACR最近基因'],['炎症/免疫相关候选',selected.filter(g=>data.all.seacr_119_Tag.near.some(a=>a[6]===g)).join('、'),selected.filter(g=>data.all.seacr_137_Tag.near.some(a=>a[6]===g)).join('、')],['全体最近基因ID数','815','4,118'],['MACS2最近基因ID数','17','72']],[280,415,415],{y:155,h:285,fontSize:23});
 text(s,'当前可汇报的功能层面结果',85,471,1110,40,28,blue,true);
 text(s,'两组基因列表均包含上述炎症相关候选，尚不能归为119特异功能。\n现有文件未提供GO/KEGG/Reactome条目、GeneRatio和校正P值。\n功能富集状态：尚未完成正式检验，不能写成“已显著富集”或“无显著富集”。',85,525,1110,133,24);
}
// 14
{
 const s=slide('结果小结',notesCommon+'仅总结本次资料中实际数据。未将染色质共占位、候选基因邻近和蛋白物理互作混同。');
 table(s,[['问题','现有观察','汇报表述'],['峰集中在哪里','MACS2中119为13/22、137为47/75触及TSS±2kb','以TSS邻近区域为主的分布'],['137与Pol II','53/75个严格共同区间相交','存在共占位候选'],['137与rRNA','少量相交，单个片段起始边缘信号较高','局部rRNA关联线索'],['119与炎症相关基因','SEACR命中若干候选，严格峰及TSS信号支持有限','存在候选基因关联，特异性未确立'],['功能富集','已有对应基因列表，缺少正式富集统计','待补GO/通路富集结果']],[250,550,310],{y:156,h:380,fontSize:23});
 text(s,'主展示：严格共同峰、同状态IgG对照和已有IGV实例。\n放松参数结果用于展示结果敏感性，蛋白互作仍需独立证据。',85,590,1110,82,25,blue,true);
}
p.delete(initial.records[0].id);
const candidate=path.join(build,'focus_candidate.pptx');
await (await PresentationFile.exportPptx(p)).save(candidate);
await fs.writeFile(path.join(build,'focus_content.json'),JSON.stringify(copy,null,2));
const preview=path.join(build,'focus_previews');await fs.mkdir(preview,{recursive:true});
for(const [i,s] of p.slides.items.entries()){const b=await p.export({slide:s,format:'png',scale:1});await fs.writeFile(path.join(preview,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await b.arrayBuffer()));}
console.log('DRAFT',candidate,'SLIDES',p.slides.items.length);
if(process.argv.includes('--draft'))process.exit(0);
const delivery=path.join(build,'deliverables');await fs.mkdir(delivery,{recursive:true});
const final=path.join(delivery,'3D4_CUTTag_区域与基因结果汇报_模板版_20261013.pptx');
const result=await finalizePresentation({workspaceDir:build,candidatePath:candidate,finalPath:final,pythonExecutable:python,
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit',...tables.flatMap(n=>['--require-native-table-slide',String(n)])],
 requiredNativeTableOwnerSlides:tables,requiredNativeChartOwnerSlides:[...new Set(charts)],materializeLiteralChartWorkbooks:true,
 fontPolicy:{basis:'design',families:[font,'Times New Roman','等线']},verifyArtifactToolImport:true,receiptPath:path.join(build,'focus_validation.json')});
await fs.mkdir(output,{recursive:true});await fs.copyFile(final,path.join(output,path.basename(final)));
console.log(JSON.stringify(result));
