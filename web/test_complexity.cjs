const test = require("node:test");
const assert = require("node:assert/strict");
const {calculate} = require("./complexity.js");
const edges=[-15,-5,5,15],xjaw=[-20,20],yjaw=[-15,15];
const cp=(l,r,m)=>({left:l,right:r,mu:m});
function near(a,b){assert.ok(Math.abs(a-b)<1e-10, String(a)+" not equal "+String(b));}
test("regular MLC has MCSv LSV AAV unity",()=>{const a=cp([-5,-5,-5],[5,5,5],0),b=cp([-5,-5,-5],[5,5,5],1);const s=calculate([a,b],edges,xjaw,yjaw);near(s.mcs,1);near(s.lsv,1);near(s.aav,1);});
test("beamAAV reduced for smaller aperture",()=>{const a=cp([-10,-10,-10],[10,10,10],0),b=cp([-5,-5,-5],[5,5,5],1);const s=calculate([a,b],edges,xjaw,yjaw);near(s.mcs,.75);near(s.lsv,1);near(s.aav,.75);});
test("jagged bank reduces LSV only",()=>{const a=cp([-5,-5,-10],[5,5,5],0),b=cp([-5,-5,-10],[5,5,5],1);const s=calculate([a,b],edges,xjaw,yjaw);near(s.mcs,.5);near(s.lsv,.5);near(s.aav,1);});
test("closed aperture returns null",()=>{const a=cp([0,0,0],[0,0,0],0),b=cp([0,0,0],[0,0,0],1);assert.equal(calculate([a,b],edges,xjaw,yjaw).mcs,null);});
test("zero MU interval not counted",()=>{const a=cp([-10,-10,-10],[10,10,10],0),b=cp([-10,-10,-10],[10,10,10],0),c=cp([-5,-5,-5],[5,5,5],1);near(calculate([a,b,c],edges,xjaw,yjaw).mcs,.75);});
test("jaw clipping affects denominator",()=>{const a=cp([-10,-10,-10],[10,10,10],0),b=cp([-1,-1,-1],[1,1,1],1);near(calculate([a,b],edges,[-3,3],[-11,11]).mcs,2/3);});
