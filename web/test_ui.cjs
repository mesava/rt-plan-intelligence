const test=require("node:test");
const assert=require("node:assert/strict");
const fs=require("node:fs");
const vm=require("node:vm");
const complexity=require("./complexity.js");
test("browser synthetic demo initializes and navigates plans/CP",()=>{
 const html=fs.readFileSync(__dirname+"/index.html","utf8");
 const scripts=[...html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)].map(x=>x[1]).filter(x=>x.trim());
 assert.equal(scripts.length,1);
 const nodes=new Map();
 const document={
   getElementById(id){
     if(!nodes.has(id))nodes.set(id,{value:id==="plan"?"pelvis":"0",max:"0",
       textContent:"",innerHTML:"",handlers:{},addEventListener(event,handler){this.handlers[event]=handler;}});
     return nodes.get(id);
   },
   querySelectorAll(){return [];}
 };
 vm.runInNewContext(scripts[0],{document,window:{RTPIComplexity:complexity}}, {timeout:3000});
 const el=id=>document.getElementById(id);
 assert.equal(el("mcp").textContent,"44");
 assert.match(el("mmcs").textContent,/0[,.][0-9]+/);
 assert.match(el("metricrows").innerHTML,/beamAAV/);
 assert.match(el("mlcSvg").innerHTML,/<rect/);
 assert.match(el("trend").innerHTML,/<polyline/);
 el("cp").value="11";el("cp").oninput();
 assert.match(el("cpbadge").textContent,/11/);
 el("plan").value="imrt";el("plan").handlers.change();
 assert.equal(el("mcp").textContent,"20");
 assert.equal(el("mmcs").textContent,"—");
 assert.equal(el("mlsv").textContent,"—");
 assert.equal(el("maav").textContent,"—");
});
