/* Synthetic-only browser CP transition diagnostics. Mirrors selected rtpi.transitions metrics.
 * Not measured MLC speed, reconstructed fluence, gamma or treatment deliverability.
 */
(function(global){
"use strict";
function calculate(cps,beamMU,sectorWidth=30){
 if(!Array.isArray(cps) || cps.length<2)
   return {transitions:[],sectors:[],summary:{maxLeafStep:0,zeroMuLeafMoves:0}};
 if(!(sectorWidth>0) || 360%sectorWidth!==0)throw Error("Invalid sector width");
 if(!Number.isFinite(beamMU)||beamMU<0)throw Error("Invalid MU");
 const final=cps[cps.length-1].mu;
 if(!(final>0))throw Error("Missing final meterset");
 const bySector=new Map(),transitions=[];
 for(let i=1;i<cps.length;i++){
  const a=cps[i-1],b=cps[i];
  if(a.left.length!==b.left.length || a.right.length!==b.right.length)
    throw Error("Mismatched leaves");
  const dc=b.mu-a.mu;
  if(!Number.isFinite(dc)||dc<0)throw Error("Decreasing CMW");
  const mu=beamMU*dc/final;
  let max=0,sum=0;
  for(let k=0;k<a.left.length;k++){
   const dl=Math.abs(b.left[k]-a.left[k]),dr=Math.abs(b.right[k]-a.right[k]);
   max=Math.max(max,dl,dr);sum+=dl+dr;
  }
  const validAngle=Number.isFinite(a.angle)&&Number.isFinite(b.angle);
  const signed=validAngle?((b.angle-a.angle+180)%360+360)%360-180:null;
  const sweep=validAngle?Math.abs(signed):null;
  const midpoint=validAngle?(a.angle+signed/2+360)%360:null;
  const row={
   from:i-1,to:i,deltaMu:mu,
   maxLeafStep:max,sumLeafStep:sum,
   maxMmPerMu:mu>0?max/mu:null,
   sweepDeg:sweep,
   maxMmPerDeg:sweep>0?max/sweep:null,
   zeroMuLeafMove:dc===0&&max>1e-8,
   angleMid:midpoint,
  };
  transitions.push(row);
  if(sweep!==null&&sweep>1e-8){
   const sector=Math.floor(midpoint/sectorWidth)*sectorWidth;
   if(!bySector.has(sector))bySector.set(sector,{start:sector,end:sector+sectorWidth,mu:0,sweepDeg:0,leafTravel:0,maxLeaf:0,count:0,zeroMuLeafMoves:0});
   const v=bySector.get(sector);
   v.mu+=mu;v.sweepDeg+=sweep;v.leafTravel+=sum;v.maxLeaf=Math.max(v.maxLeaf,max);v.count++;
   if(row.zeroMuLeafMove)v.zeroMuLeafMoves++;
  }
 }
 const sectors=[...bySector.values()].sort((a,b)=>a.start-b.start).map(v=>({...v,muPerDegree:v.sweepDeg?v.mu/v.sweepDeg:null,leafStepPerMu:v.mu?v.leafTravel/v.mu:null}));
 return {
  transitions,sectors,
  summary:{
   maxLeafStep:Math.max(0,...transitions.map(t=>t.maxLeafStep)),
   maxMmPerMu:Math.max(0,...transitions.map(t=>t.maxMmPerMu??0)),
   zeroMuLeafMoves:transitions.filter(t=>t.zeroMuLeafMove).length,
   sectorCount:sectors.length,
  },
 };
}
const api={calculate};
if(typeof module!=="undefined" && module.exports)module.exports=api;
global.RTPITransitions=api;
})(typeof window!=="undefined"?window:globalThis);
