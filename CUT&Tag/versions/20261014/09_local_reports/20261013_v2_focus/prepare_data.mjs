import fs from 'node:fs/promises';
import path from 'node:path';
const root='C:/Users/Luca/Documents/Codex/2026-09-17/cut-tag';
const screen=path.join(root,'work/hypothesis_screen_import/hypothesis_screen_current_peaks');
const signal=path.join(root,'work/hypothesis_signal_qc_import_51892/job-51892');
const build=path.join(root,'work/ppt_focus_20261013');
const raw=async p=>(await fs.readFile(p,'utf8')).trim().split(/\r?\n/).filter(Boolean).map(x=>x.split('\t'));
const objects=async p=>{const [h,...r]=await raw(p);return r.map(a=>Object.fromEntries(h.map((k,i)=>[k,a[i]])));};
const bed=async p=>(await raw(p)).map(a=>({chr:a[0],start:+a[1],end:+a[2],fields:a}));
const genes=await bed(path.join(screen,'annotations/Ensembl115.genes.bed'));
const bychr=new Map();for(const g of genes){if(!bychr.has(g.chr))bychr.set(g.chr,[]);bychr.get(g.chr).push(g);}
const overlap=(a,b)=>a.chr===b.chr&&a.start<b.end&&b.start<a.end;
const all={};
for(const caller of ['macs2','seacr'])for(const state of ['119','137'])for(const mark of ['Tag','Pol2']){
 const key=`${caller}_${state}_${mark}`;
 const peaks=await bed(path.join(screen,`peak_sets/${key}_rep1-rep2_shared.bed`));
 const near=await raw(path.join(screen,`nearest_genes/${key}.nearest_gene.tsv`));
 const nearGenes=[...new Set(near.filter(a=>a[7]!=='.'&&a[7]!=='-1').map(a=>a[7]))];
 const regions={promoter:0,body:0,distal:0};
 const chr={};for(const p of peaks){chr[p.chr]=(chr[p.chr]||0)+1;let prom=false,body=false;
  for(const g of bychr.get(p.chr)||[]){const t=g.fields[5]==='-'?g.end-1:g.start;
   if(p.start<t+2000&&t-2000<p.end)prom=true;
   if(overlap(p,g))body=true;
  }
  regions[prom?'promoter':body?'body':'distal']++;
 }
 all[key]={peaks,near,nearGeneCount:nearGenes.length,regions,chr};
}
const compare={};for(const c of ['macs2','seacr']){
 for(const st of ['119','137']){const a=all[`${c}_${st}_Tag`].peaks,b=all[`${c}_${st}_Pol2`].peaks;compare[`${c}_${st}_Pol2`]={total:a.length,hit:a.filter(p=>b.some(q=>overlap(p,q))).length};}
 const a=all[`${c}_119_Tag`].peaks,b=all[`${c}_137_Tag`].peaks;
 compare[`${c}_119_vs137`]={total:a.length,hit:a.filter(p=>b.some(q=>overlap(p,q))).length};
}
const tss=await objects(path.join(signal,'candidate_TSS_plusminus500bp.exact_signal.tsv'));
const rrna=await objects(path.join(signal,'rRNA.exact_signal.tsv'));
const rrnaBins=(await objects(path.join(signal,'rRNA_and_flanks_100bp_bins.exact_signal.tsv'))).filter(r=>r.chrom==='AEMK02000489.1'&&r.name.includes('43231-45101'));
const candidates={};for(const c of ['macs2','seacr'])for(const st of ['119','137']){
 const p=path.join(screen,`overlaps/${c}_${st}_Tag.candidate_gene_window_hits.tsv`);
 const str=(await fs.readFile(p,'utf8')).trim();candidates[`${c}_${st}`]=str?[...new Set(str.split(/\r?\n/).map(r=>r.split('\t')[3]))]:[];
}
const data={all,compare,tss,rrna,rrnaBins,candidates};
await fs.writeFile(path.join(build,'data.json'),JSON.stringify(data));
console.log(JSON.stringify({sets:Object.fromEntries(Object.entries(all).filter(([k])=>k.endsWith('Tag')).map(([k,v])=>[k,{n:v.peaks.length,genes:v.nearGeneCount,regions:v.regions,chrom:v.chr}])),compare,candidates,rRNAbins:rrnaBins.length},null,2));
