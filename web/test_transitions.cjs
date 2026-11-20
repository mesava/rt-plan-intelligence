const test=require('node:test');
const assert=require('node:assert/strict');
const {calculate}=require('./transitions.js');
const cp=(left,right,mu,angle)=>({left,right,mu,angle});
test('shortest arc and normalized MU',()=>{
 const a=cp([-5,-5],[5,5],0,359),b=cp([-3,-3],[7,7],.5,1),c=cp([-2,-2],[8,8],1,3);
 const r=calculate([a,b,c],120,30);
 assert.equal(r.transitions.length,2);
 assert.equal(r.transitions[0].deltaMu,60);
 assert.equal(r.transitions[0].sweepDeg,2);
 assert.equal(r.transitions[0].maxLeafStep,2);
 assert.equal(r.sectors.length,1);
 assert.equal(r.sectors[0].start,0);
 assert.equal(r.sectors[0].mu,120);
});
test('stationary IMRT has no angular sectors',()=>{
 const a=cp([-5,-5],[5,5],0,100),b=cp([-4,-4],[5,5],1,100);
 const r=calculate([a,b],100);
 assert.equal(r.sectors.length,0);
 assert.equal(r.transitions[0].maxMmPerDeg,null);
 assert.equal(r.transitions[0].maxMmPerMu,.01);
});
test('zero MU transition flagged separately',()=>{
 const a=cp([-5],[5],0,10),b=cp([-3],[7],0,11),c=cp([-3],[7],1,12);
 const r=calculate([a,b,c],100);
 assert.equal(r.transitions[0].zeroMuLeafMove,true);
 assert.equal(r.transitions[0].maxMmPerMu,null);
 assert.equal(r.summary.zeroMuLeafMoves,1);
});
test('cannot treat negative CMW as valid',()=>{
 const a=cp([-5],[5],1,10),b=cp([-5],[5],0,12);
 assert.throws(()=>calculate([a,b],100),/Decreasing CMW/);
});
