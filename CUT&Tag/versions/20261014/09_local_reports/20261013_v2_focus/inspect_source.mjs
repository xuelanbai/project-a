import fs from 'node:fs/promises';
import path from 'node:path';
import {FileBlob, PresentationFile} from 'file:///C:/Users/Luca/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';
const root='C:/Users/Luca/Documents/Codex/2026-09-17/cut-tag/work/ppt_focus_20261013';
await fs.mkdir(root,{recursive:true});
for(const [label,file] of [
 ['reference','D:/Lucas临时文件/会议汇报/工作总结-10.pptx'],
 ['previous','C:/Users/Luca/Documents/Codex/2026-09-17/cut-tag/outputs/ppt_report_20261013/3D4_CUTTag_结果汇报_简洁版_20261013.pptx']
]){
 const p=await PresentationFile.importPptx(await FileBlob.load(file));
 const dir=path.join(root,label);await fs.mkdir(dir,{recursive:true});
 const inspect=await p.inspect({select:'$.slides[*]',include:['id','index','layout','title','elements.id','elements.kind','elements.text','elements.position','elements.style','notes.text']});
 await fs.writeFile(path.join(dir,'inspect.ndjson'),inspect.ndjson);
 console.log(label,'slides',p.slides.items.length,'masters',p.masters.items.length,'layouts',p.layouts.items.length);
 for(const [i,s] of p.slides.items.entries()){
  const png=await p.export({slide:s,format:'png',scale:1});
  await fs.writeFile(path.join(dir,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await png.arrayBuffer()));
 }
 console.log('rendered',label);
}
