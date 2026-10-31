/* Research-only, synthetic MCSv / LSV / beamAAV browser implementation. */
(function (global) {
"use strict";
function calculate(cps,edges,xjaw,yjaw){
  if(!Array.isArray(cps) || cps.length === 0) return {mcs:null,lsv:null,aav:null,lsvByCP:[],aavByCP:[]};
  const n=edges.length-1;
  if(n<1 || cps.some(c=>c.left.length!==n || c.right.length!==n)) throw new Error("Invalid leaf geometry");
  const heights=edges.slice(0,-1).map((lo,i)=>Math.max(0,Math.min(edges[i+1],yjaw[1])-Math.max(lo,yjaw[0])));
  const ll=cps.map(c=>c.left.map(x=>Math.max(x,xjaw[0])));
  const rr=cps.map(c=>c.right.map(x=>Math.min(x,xjaw[1])));
  let envelope=0;
  for(let i=0;i<n;i++) envelope+=Math.max(0,Math.max(...rr.map(p=>p[i]))-Math.min(...ll.map(p=>p[i])))*heights[i];
  function bankScore(values){
    if(values.length<2)return 1;
    const spread=Math.max(...values)-Math.min(...values);
    if(spread===0)return 1;
    let movement=0;
    for(let i=1;i<values.length;i++)movement+=Math.abs(values[i]-values[i-1]);
    return Math.max(0,Math.min(1,1-movement/((values.length-1)*spread)));
  }
  const lsvByCP=[],aavByCP=[];
  for(const cp of cps){
    let area=0;const L=[],R=[];
    for(let i=0;i<n;i++){
      const l=Math.max(cp.left[i],xjaw[0]),r=Math.min(cp.right[i],xjaw[1]);
      const gap=Math.max(0,r-l);area+=gap*heights[i];
      if(gap>0 && heights[i]>0){L.push(cp.left[i]);R.push(cp.right[i]);}
    }
    lsvByCP.push(bankScore(L)*bankScore(R));
    aavByCP.push(envelope>0?area/envelope:null);
  }
  if(cps.length<2 || !(envelope>0))return {mcs:null,lsv:null,aav:null,lsvByCP,aavByCP};
  let total=0,mcs=0,lsv=0,aav=0;
  for(let i=1;i<cps.length;i++){
    const dm=cps[i].mu-cps[i-1].mu;
    if(!Number.isFinite(dm) || dm<0)throw new Error("Invalid meterset");
    if(dm<=0)continue;
    const li=(lsvByCP[i-1]+lsvByCP[i])/2;
    const ai=(aavByCP[i-1]+aavByCP[i])/2;
    total+=dm;mcs+=dm*li*ai;lsv+=dm*li;aav+=dm*ai;
  }
  if(!(total>0))return {mcs:null,lsv:null,aav:null,lsvByCP,aavByCP};
  return {mcs:mcs/total,lsv:lsv/total,aav:aav/total,lsvByCP,aavByCP};
}
const api={calculate};
if(typeof module!=="undefined" && module.exports)module.exports=api;
global.RTPIComplexity=api;
})(typeof window!=="undefined"?window:globalThis);
