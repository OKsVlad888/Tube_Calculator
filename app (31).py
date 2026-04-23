"""
High-Pressure Gas Flow Calculator — Full UI Edition v2
=======================================================
BUG FIXES vs previous version:
  1. HUD embedded as base64 inline — no hud_background.html needed on server
  2. toKey() fixed: parens removed (not replaced with _) — fields now appear
  3. Math.pow() replaces ** operator — full browser compatibility
  4. getCurrentInletP uses toKey() — correct live theme updates
  5. Safety block (red) / warning (orange) fully wired
  6. Confidence bar displayed after every successful calculation

Run:  streamlit run app.py   (single file, no external dependencies)
"""

import base64
import json
import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(
    page_title="High-Pressure Gas Flow Calculator",
    page_icon="💨",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown("""
<style>
  #MainMenu,footer,header,
  [data-testid="stToolbar"],[data-testid="stDecoration"],
  [data-testid="stStatusWidget"],[data-testid="collapsedControl"]{display:none!important}
  .stApp{background:#000d0d!important}
  section.main,.block-container{padding:0!important;margin:0!important;max-width:100vw!important}
  iframe[title="streamlit_components_v1.html_v1"]{
    position:fixed!important;inset:0!important;
    width:100vw!important;height:100vh!important;
    border:none!important;z-index:100!important;
    display:block!important}
</style>
""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════════════
#  TUBE-SPEC TABLE — DO NOT MODIFY DATA OR LOGIC
# ══════════════════════════════════════════════════════════════════════════════
TUBE_TABLE = [
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":4.9,    "tube_od_w":'1/4"X0.028"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":7.036,  "tube_od_w":'3/8"X0.049"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":9.398,  "tube_od_w":'1/2"X0.065"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":12.573, "tube_od_w":'5/8"X0.065"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":14.834, "tube_od_w":'3/4"X0.083"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":17.399, "tube_od_w":'7/8"X0.095"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":19.863, "tube_od_w":'1"X0.109"',  "max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":3.86,   "tube_od_w":'1/4"X0.049"',    "max_pressure":400.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":6.223,  "tube_od_w":'3/8"X0.065"',    "max_pressure":400.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":8.468,  "tube_od_w":'1/2"X0.083"',    "max_pressure":400.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":9.119,  "tube_od_w":'9/16"X0.10175"', "max_pressure":400.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":13.106, "tube_od_w":'3/4"X0.117"',    "max_pressure":400.0},
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S11","id_mm":2.108,  "tube_od_w":'1/4"X0.083"',"max_pressure":2000.0},
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S11","id_mm":3.1753, "tube_od_w":'3/8"X0.125"',"max_pressure":2000.0},
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S12","id_mm":2.6786, "tube_od_w":'1/4"X0.705"',   "max_pressure":1380.0},
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S12","id_mm":5.516,  "tube_od_w":'3/8"X0.0860"',  "max_pressure":1380.0},
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S12","id_mm":7.925,  "tube_od_w":'9/16"X0.12525"',"max_pressure":1380.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":4.572,  "tube_od_w":'1/4"X0.035"',"max_pressure":238.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":7.747,  "tube_od_w":'3/8"X0.035"',"max_pressure":170.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":10.21,  "tube_od_w":'1/2"X0.049"',"max_pressure":170.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":15.748, "tube_od_w":'3/4"X0.065"',"max_pressure":170.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":22.1,   "tube_od_w":'1"X0.065"',  "max_pressure":102.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":34.8,   "tube_od_w":'1.5"X0.065"',"max_pressure":68.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":47.5,   "tube_od_w":'2"X0.065"',  "max_pressure":61.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":4.572,    "tube_od_w":'1/4"X0.035"', "max_pressure":204.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":7.747,    "tube_od_w":'3/8"X0.035"', "max_pressure":170.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":10.21,    "tube_od_w":'1/2"X0.049"', "max_pressure":170.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":13.38585, "tube_od_w":'5/8"X0.049"', "max_pressure":170.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":15.748,   "tube_od_w":'3/4"X0.065"', "max_pressure":170.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":4.572,   "tube_od_w":'1/4"X0.035"', "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":7.0358,  "tube_od_w":'3/8"X0.049"', "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":10.2108, "tube_od_w":'1/2"X0.049"', "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":13.3858, "tube_od_w":'5/8"X0.049"', "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":15.748,  "tube_od_w":'3/4"X0.065"', "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":18.923,  "tube_od_w":'7/8"X0.065"', "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":22.098,  "tube_od_w":'1"X0.065"',   "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":28.448,  "tube_od_w":'1.25"X0.065"',"max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":34.4424, "tube_od_w":'1.5"X0.072"', "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":46.5836, "tube_od_w":'2"X0.083"',   "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":58.674,  "tube_od_w":'2.5"X0.095"', "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":70.6628, "tube_od_w":'3"X0.109"',   "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"B1","id_mm":94.7928, "tube_od_w":'4"X0.134"',   "max_pressure":40.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":11.6, "tube_od_w":"16mmX2.2mm",  "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":14.4, "tube_od_w":"20mmX2.8mm",  "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":18.0, "tube_od_w":"25mmX3.5mm",  "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":23.2, "tube_od_w":"32mmX4.4mm",  "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":29.0, "tube_od_w":"40mmX5.5mm",  "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":36.2, "tube_od_w":"50mmX6.9mm",  "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":45.8, "tube_od_w":"63mmX8.6mm",  "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":54.4, "tube_od_w":"75mmX10.3mm", "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":65.4, "tube_od_w":"90mmX12.3mm", "max_pressure":10.0},
    {"gas_codes":["N2","Ar","Air"],"spec":"P6","id_mm":79.8, "tube_od_w":"110mmX15.1mm","max_pressure":10.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":15.8,  "tube_od_w":'1/2" SCH40', "max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":15.59, "tube_od_w":'3/4" SCH40', "max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":26.24, "tube_od_w":'1" SCH40',   "max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":35.04, "tube_od_w":'1.25" SCH40',"max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":40.9,  "tube_od_w":'1.5" SCH40', "max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":52.51, "tube_od_w":'2" SCH40',   "max_pressure":15.0},
]

TUBE_TABLE_JSON = json.dumps(TUBE_TABLE, ensure_ascii=False)

# ══════════════════════════════════════════════════════════════════════════════
#  HUD BACKGROUND — EMBEDDED INLINE as base64 (no hud_background.html needed)
#  The entire animated HUD is encoded here so Streamlit Cloud doesn't need
#  any extra files. If hud_background.html exists in the same directory,
#  it is used instead (lets you iterate on the HUD without redeploying).
# ══════════════════════════════════════════════════════════════════════════════
_HUD_HTML = b"""<!DOCTYPE html>
<html><head><meta charset="UTF-8">
<style>
*{margin:0;padding:0;box-sizing:border-box}
html,body{width:100%;height:100%;overflow:hidden;background:#000d0d;font-family:'Courier New',monospace}
.gb{position:fixed;inset:0;background-image:linear-gradient(rgba(0,229,204,.03) 1px,transparent 1px),linear-gradient(90deg,rgba(0,229,204,.03) 1px,transparent 1px);background-size:50px 50px;pointer-events:none;z-index:0}
.sl{position:fixed;inset:0;background:repeating-linear-gradient(to bottom,transparent 0,transparent 3px,rgba(0,229,204,.01) 3px,rgba(0,229,204,.01) 4px);pointer-events:none;z-index:1}
#hc{position:fixed;inset:0;z-index:2;pointer-events:none}
.co{position:fixed;width:50px;height:50px;border-color:rgba(0,229,204,.5);border-style:solid;z-index:3;pointer-events:none}
.tl{top:8px;left:8px;border-width:2px 0 0 2px}.tr{top:8px;right:8px;border-width:2px 2px 0 0}
.bl{bottom:8px;left:8px;border-width:0 0 2px 2px}.br{bottom:8px;right:8px;border-width:0 2px 2px 0}
.hu{position:fixed;font-family:'Courier New',monospace;color:rgba(0,229,204,.7);pointer-events:none;z-index:3}
.hl{font-size:8px;letter-spacing:2px;text-transform:uppercase;color:rgba(0,229,204,.38);display:block}
.hv{font-size:20px;font-weight:bold;color:#00ffee;display:block;text-shadow:0 0 8px rgba(0,229,204,.5)}
.sp{position:fixed;left:14px;top:50%;transform:translateY(-50%);font-size:32px;color:rgba(0,229,204,.28);animation:sp 12s linear infinite;z-index:3;pointer-events:none}
@keyframes sp{from{transform:translateY(-50%) rotate(0)}to{transform:translateY(-50%) rotate(360deg)}}
.bl2{animation:bl 0.9s step-end infinite}@keyframes bl{50%{opacity:0}}
.xbt{display:inline-block;width:9px;height:48px;background:rgba(0,229,204,.07);border:1px solid rgba(0,229,204,.16);vertical-align:bottom;margin:0 1px;position:relative;overflow:hidden}
.xbf{position:absolute;bottom:0;left:0;width:100%;background:rgba(0,229,204,.62)}
</style></head><body>
<div class="gb"></div><div class="sl"></div>
<div class="co tl"></div><div class="co tr"></div><div class="co bl"></div><div class="co br"></div>
<div class="sp">&#10022;</div>
<div class="hu" style="top:65px;left:18px;font-size:10px">
  <span class="hl">Data Channel Identification</span><span class="hv" id="v1">89</span>
  <span class="bl2" style="color:#00e5cc;font-size:10px">&#11044;</span>
  <div style="margin-top:6px"><span class="hl">Data Channel Identification</span><span class="hv" id="v2">33</span></div>
  <div style="display:flex;gap:10px;font-size:11px;color:#00e5cc;margin-top:2px">
    <span id="v3">31</span><span style="opacity:.55">29</span><span style="opacity:.55">29</span>
  </div>
  <div style="display:flex;gap:10px;font-size:10px;color:rgba(0,229,204,.45)"><span>26</span><span>57</span><span>39</span></div>
</div>
<div class="hu" style="top:158px;left:18px;font-size:8px">
  <span class="hl">Access Point</span>
  <div style="font-size:8px;color:rgba(0,229,204,.38);max-width:110px;margin-top:2px">Lorem ipsum dolor sit amet lorem</div>
  <div style="margin-top:4px;font-size:8px;color:rgba(0,229,204,.28)">SCANNING DATA <span class="bl2">&#9646;</span><br>SCANNING DATA</div>
</div>
<div class="hu" style="top:10px;left:44%;font-size:9px;white-space:nowrap">
  <span class="hl">Data Channel Identification</span>
  <div style="display:flex;gap:3px;align-items:flex-end;margin-top:2px">
    <span style="font-size:12px;color:#00ffee;margin-right:4px" id="vtc">26</span>
    <div class="xbt"><div class="xbf" id="xb0" style="height:60%"></div></div>
    <div class="xbt"><div class="xbf" id="xb1" style="height:40%"></div></div>
    <div class="xbt"><div class="xbf" id="xb2" style="height:75%"></div></div>
    <div class="xbt"><div class="xbf" id="xb3" style="height:50%"></div></div>
    <div class="xbt"><div class="xbf" id="xb4" style="height:85%"></div></div>
    <div class="xbt"><div class="xbf" id="xb5" style="height:35%"></div></div>
    <div class="xbt"><div class="xbf" id="xb6" style="height:65%"></div></div>
    <span style="font-size:9px;color:rgba(0,229,204,.48);margin-left:4px">&#9608; EXTRACT DATA</span>
  </div>
  <div style="display:flex;gap:15px;font-size:10px;color:rgba(0,229,204,.45);margin-top:2px"><span>29</span><span>57</span><span>39</span></div>
</div>
<div class="hu" style="top:65px;right:18px;text-align:right;max-width:190px">
  <div style="font-size:17px;font-weight:bold;color:#00ffee;letter-spacing:3px;text-shadow:0 0 9px rgba(0,229,204,.5)">LOREM IPSUM</div>
  <div style="font-size:8px;color:rgba(0,229,204,.35);margin-top:2px;max-width:185px;text-align:right">Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod</div>
  <div style="margin-top:3px;font-size:9px;color:rgba(0,229,204,.55)"><span style="color:#00e5cc">&#9658;&#9658;&#9658;</span> DATA X54</div>
</div>
<div class="hu" style="top:50%;right:18px;transform:translateY(-50%);text-align:right">
  <span class="hv" id="vm1" style="font-size:28px">75</span>
  <div style="font-size:7px;color:rgba(0,229,204,.35);letter-spacing:1.5px;text-transform:uppercase;margin-top:1px">Data Channel<br>Identification</div>
  <div style="margin:6px 0;display:flex;justify-content:flex-end">
    <svg width="44" height="44" viewBox="0 0 44 44">
      <circle cx="22" cy="22" r="18" fill="none" stroke="rgba(0,229,204,.11)" stroke-width="4"/>
      <circle cx="22" cy="22" r="18" fill="none" stroke="#00e5cc" stroke-width="4"
        stroke-dasharray="56 56" stroke-dashoffset="14" stroke-linecap="round" id="rarc"
        style="filter:drop-shadow(0 0 4px rgba(0,229,204,.45))"/>
    </svg>
  </div>
  <span class="hv" id="vm2" style="font-size:28px">48</span>
  <div style="font-size:7px;color:rgba(0,229,204,.35);letter-spacing:1.5px;text-transform:uppercase;margin-top:1px">Data Channel<br>Identification</div>
</div>
<div class="hu" style="bottom:28px;left:18px;font-size:8px;color:rgba(0,229,204,.35);width:175px;line-height:1.4">
  <span class="hl" style="margin-bottom:2px">Tracking Target 07DC9EA1</span>
  <div id="bs"></div>
</div>
<canvas id="hc"></canvas>
<script>
var cvs=document.getElementById("hc"),ctx=cvs.getContext("2d"),W,H,t=0;
function rsz(){W=cvs.width=window.innerWidth;H=cvs.height=window.innerHeight}
rsz();window.addEventListener("resize",rsz);
var TC="rgba(0,229,204,",TC2="#00e5cc";
var bL=[],bD=document.getElementById("bs");
var HX="0123456789ABCDEFabcdef ";
function rL(){var s="";for(var i=0;i<17;i++)s+=HX[Math.floor(Math.random()*HX.length)];return s}
for(var i=0;i<10;i++)bL.push(rL());
setInterval(function(){bL.shift();bL.push(rL());bD.textContent=bL.join("\n")},140);
var ct2=[
  {el:document.getElementById("v1"),v:89,t:89},
  {el:document.getElementById("v2"),v:33,t:33},
  {el:document.getElementById("v3"),v:31,t:31},
  {el:document.getElementById("vm1"),v:75,t:75},
  {el:document.getElementById("vm2"),v:48,t:48},
  {el:document.getElementById("vtc"),v:26,t:26}
];
setInterval(function(){ct2.forEach(function(c){if(Math.random()<0.07)c.t=10+Math.floor(Math.random()*90);c.v+=(c.t-c.v)*0.1;if(c.el)c.el.textContent=Math.round(c.v)})},200);
var xH=[60,40,75,50,85,35,65],xT=xH.slice();
setInterval(function(){for(var i=0;i<7;i++){if(Math.random()<0.25)xT[i]=8+Math.random()*92;xH[i]+=(xT[i]-xH[i])*0.1;var e=document.getElementById("xb"+i);if(e)e.style.height=xH[i].toFixed(1)+"%"}},90);
var ro=14;
setInterval(function(){ro=(ro+0.5)%113;var a=document.getElementById("rarc");if(a)a.setAttribute("stroke-dashoffset",ro.toFixed(1))},30);
var bars=[],bT=[];for(var i=0;i<9;i++){var v=0.2+Math.random()*0.8;bars.push(v);bT.push(v)}
setInterval(function(){for(var i=0;i<bars.length;i++){if(Math.random()<0.2)bT[i]=0.08+Math.random()*0.92;bars[i]+=(bT[i]-bars[i])*0.05}},70);
var sb=[0.5,0.7,0.4,0.8,0.6,0.9,0.3,0.75],sbT=sb.slice();
setInterval(function(){for(var i=0;i<sb.length;i++){if(Math.random()<0.15)sbT[i]=0.1+Math.random()*0.9;sb[i]+=(sbT[i]-sb[i])*0.05}},80);
var hl=[];
function gh(){var s="",hx="0123456789ABCDEF";for(var c=0;c<7;c++){s+=hx[Math.floor(Math.random()*16)]+hx[Math.floor(Math.random()*16)];if(c<6)s+=" "}return s}
for(var i=0;i<10;i++)hl.push(gh());
setInterval(function(){hl[Math.floor(Math.random()*hl.length)]=gh()},220);
var dots=[];for(var i=0;i<26;i++)dots.push({bx:.72+Math.random()*.22,by:.50+Math.random()*.24,ph:Math.random()*6.28,r:1.5+Math.random()*2.5,sp:.008+Math.random()*.012});
function lb(x,y,s,sz){ctx.font=(sz||8)+'px "Courier New"';ctx.fillStyle=TC+".36)";ctx.fillText(s.toUpperCase(),x,y)}
function dw(rx,ry,rw,rh,fr,am,ph,sp,col,lw){
  var x=rx*W,y=ry*H,ww=rw*W,hh=rh*H;
  ctx.fillStyle="rgba(0,18,16,.30)";ctx.fillRect(x,y,ww,hh);ctx.strokeStyle=TC+".09)";ctx.lineWidth=1;ctx.strokeRect(x,y,ww,hh);
  for(var i=1;i<4;i++){ctx.beginPath();ctx.strokeStyle=TC+".06)";ctx.moveTo(x+i*ww/4,y);ctx.lineTo(x+i*ww/4,y+hh);ctx.stroke()}
  ctx.beginPath();ctx.strokeStyle=col;ctx.lineWidth=lw;
  for(var i=0;i<=230;i++){var px=x+(i/230)*ww,vl=0;for(var f=0;f<fr.length;f++)vl+=Math.sin(i/230*Math.PI*2*fr[f]+ph[f]+t*sp[f]*60)*am[f];var py=y+hh/2+vl*(hh/2)*.82;i===0?ctx.moveTo(px,py):ctx.lineTo(px,py)}
  ctx.shadowBlur=7;ctx.shadowColor=TC2;ctx.stroke();ctx.shadowBlur=0}
function dm(){
  var x=.575*W,y=.34*H,ww=.385*W,hh=.15*H;
  ctx.fillStyle="rgba(0,18,16,.26)";ctx.fillRect(x,y,ww,hh);ctx.strokeStyle=TC+".08)";ctx.lineWidth=1;ctx.strokeRect(x,y,ww,hh);
  for(var i=1;i<4;i++){ctx.beginPath();ctx.strokeStyle=TC+".05)";ctx.moveTo(x+i*ww/4,y);ctx.lineTo(x+i*ww/4,y+hh);ctx.stroke();ctx.beginPath();ctx.moveTo(x,y+i*hh/4);ctx.lineTo(x+ww,y+i*hh/4);ctx.stroke()}
  ctx.beginPath();ctx.strokeStyle=TC+".66)";ctx.lineWidth=1.2;
  for(var i=0;i<=260;i++){var px=x+(i/260)*ww,py=y+hh/2+Math.sin(i/260*Math.PI*3.5+t*.04*60)*hh*.22+Math.sin(i/260*Math.PI*7+t*.065*60+.8)*hh*.09+Math.sin(i/260*Math.PI*1.3+t*.018*60+2)*hh*.05;i===0?ctx.moveTo(px,py):ctx.lineTo(px,py)}
  ctx.shadowBlur=5;ctx.shadowColor=TC2;ctx.stroke();ctx.shadowBlur=0;lb(x+4,y+10,"Scanning Data",7)}
function ds(){dots.forEach(function(d){var dx=d.bx*W+Math.sin(d.ph+t*d.sp*60)*6,dy=d.by*H+Math.cos(d.ph*1.4+t*d.sp*50)*5,a=.28+.45*Math.sin(d.ph+t*.03*60);ctx.beginPath();ctx.arc(dx,dy,d.r,0,6.28);ctx.fillStyle=TC+a+")";ctx.shadowBlur=5;ctx.shadowColor=TC2;ctx.fill();ctx.shadowBlur=0})}
function db(){
  var bx=.735*W,by=.78*H,bw=.23*W,bh=.155*H;
  ctx.fillStyle="rgba(0,18,16,.24)";ctx.fillRect(bx,by,bw,bh);ctx.strokeStyle=TC+".08)";ctx.lineWidth=1;ctx.strokeRect(bx,by,bw,bh);lb(bx+3,by-4,"Analysed OK",7);
  var slot=bw/bars.length;bars.forEach(function(v,i){var bX=bx+i*slot+slot*.15,bW=slot*.7,bH=v*bh*.86,bY=by+bh-bH;var g=ctx.createLinearGradient(0,bY,0,by+bh);g.addColorStop(0,TC+(0.32+v*.5)+")");g.addColorStop(1,TC+".06)");ctx.fillStyle=g;ctx.fillRect(bX,bY,bW,bH)})}
function dsb(){
  var bx=.875*W,by=.73*H,bw=.108*W,bh=.10*H;
  ctx.fillStyle="rgba(0,18,16,.23)";ctx.fillRect(bx,by,bw,bh);ctx.strokeStyle=TC+".08)";ctx.lineWidth=1;ctx.strokeRect(bx,by,bw,bh);lb(bx+2,by-3,"Get Data",7);
  var slot=bw/sb.length;sb.forEach(function(v,i){var bX=bx+i*slot+slot*.12,bW=slot*.76,bH=v*bh*.86,bY=by+bh-bH;ctx.fillStyle=TC+(0.26+v*.5)+")";ctx.fillRect(bX,bY,bW,bH)})}
function dhex(){var x=.785*W,y=.44*H;ctx.font='7.5px "Courier New"';ctx.fillStyle=TC+".26)";hl.forEach(function(l,i){ctx.fillText(l,x,y+i*10.5)})}
function dvig(){var g=ctx.createRadialGradient(W/2,H/2,H*.22,W/2,H/2,H*.82);g.addColorStop(0,"transparent");g.addColorStop(1,"rgba(0,0,0,.5)");ctx.fillStyle=g;ctx.fillRect(0,0,W,H)}
function draw(){
  ctx.clearRect(0,0,W,H);dvig();
  ctx.strokeStyle=TC+".11)";ctx.lineWidth=1;
  ctx.beginPath();ctx.moveTo(65,13);ctx.lineTo(W-65,13);ctx.stroke();
  ctx.beginPath();ctx.moveTo(65,H-13);ctx.lineTo(W-65,H-13);ctx.stroke();
  dw(.55,.06,.42,.10,[2.1,4.5],[.35,.14],[0,1.2],[.028,.05],TC+".70)",1.4);
  dw(.55,.19,.42,.09,[3.2,6.0],[.30,.11],[.5,2.0],[.022,.06],TC+".54)",1.2);
  dm();ds();db();dsb();dhex();
  ctx.font='7.5px "Courier New"';ctx.fillStyle=TC+".30)";
  ctx.fillText("ANALYSING THE",.60*W,.30*H);ctx.fillText("SCANNING DATA",.60*W,.30*H+10);
  ctx.fillText("ANALYSED DATA",.60*W,.76*H);
  ctx.fillText("IDENTIFICATION",.86*W,.30*H);ctx.fillText("READY",.86*W,.30*H+10);
  t+=1/60;requestAnimationFrame(draw)}
draw();
</script></body></html>"""

# Build the HUD data URI inline — works on Streamlit Cloud with no extra files
_HUD_B64 = "data:text/html;base64," + base64.b64encode(_HUD_HTML).decode()
HUD_TAG  = f"<iframe id='hud' src='{_HUD_B64}' sandbox='allow-scripts'></iframe>"

# ══════════════════════════════════════════════════════════════════════════════
#  FULL PAGE HTML
# ══════════════════════════════════════════════════════════════════════════════
PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
html,body{width:100%;height:100%;overflow:hidden;background:#000d0d;
  font-family:'Courier New',monospace;font-size:14px;color:#00e5cc}

/* HUD sits behind everything */
#hud{position:fixed;inset:0;width:100%;height:100%;border:none;pointer-events:none;z-index:1}

/* Outer wrapper */
#wrapper{position:fixed;inset:0;display:flex;align-items:center;justify-content:center;
  z-index:500;pointer-events:none;padding:10px}

/* Main panel — NO internal scroll */
#panel{
  pointer-events:all;width:min(98%,552px);max-height:97vh;overflow:hidden;
  background:var(--pbg);border:1.5px solid var(--tc);border-radius:4px;
  padding:10px 20px 14px;position:relative;
  animation:pglow 2.5s ease-in-out infinite alternate;
  transition:background .35s,border-color .35s}
@keyframes pglow{
  from{box-shadow:0 0 8px var(--gs),inset 0 0 20px rgba(0,229,204,.04)}
  to  {box-shadow:0 0 30px var(--ge),inset 0 0 45px rgba(0,229,204,.07)}}
#panel::after{content:"";position:absolute;inset:0;pointer-events:none;border-radius:4px;
  background:repeating-linear-gradient(to bottom,transparent 0,transparent 3px,
    rgba(0,229,204,.007) 3px,rgba(0,229,204,.007) 4px)}

/* ── CSS theme variables ────────────────────────────────────────────────── */
:root{
  --tc:#00e5cc;--tc2:#00ffee;--tc3:rgba(0,229,204,.36);
  --pbg:rgba(0,8,8,.95);--ibg:rgba(0,14,12,.98);--lc:rgba(0,229,204,.70);
  --bbg:rgba(0,229,204,.07);--bh:rgba(0,229,204,.20);
  --gs:rgba(0,229,204,.25);--ge:rgba(0,229,204,.66);--ashw:none}
/* Orange warning theme */
.W{--tc:#ffb800;--tc2:#ffe590;--tc3:rgba(255,184,0,.33);
   --pbg:rgba(18,10,0,.97);--ibg:rgba(22,14,0,.99);--lc:rgba(255,200,80,.72);
   --bbg:rgba(255,184,0,.07);--bh:rgba(255,184,0,.22);
   --gs:rgba(255,184,0,.25);--ge:rgba(255,184,0,.68);--ashw:flex}
/* Red danger theme */
.D{--tc:#ff3333;--tc2:#ff8888;--tc3:rgba(255,60,60,.33);
   --pbg:rgba(18,0,0,.98);--ibg:rgba(22,2,2,.99);--lc:rgba(255,140,140,.72);
   --bbg:rgba(255,60,60,.07);--bh:rgba(255,60,60,.22);
   --gs:rgba(255,60,60,.25);--ge:rgba(255,60,60,.68);--ashw:flex}

/* ── Alert header ──────────────────────────────────────────────────────── */
#alhdr{display:var(--ashw);align-items:center;justify-content:center;gap:8px;
  padding:5px 10px;margin-bottom:7px;background:rgba(0,0,0,.18);
  border:1px solid var(--tc3);border-radius:2px;position:relative;z-index:1;
  animation:apls 1.3s ease-in-out infinite alternate;transition:border-color .35s}
@keyframes apls{from{box-shadow:0 0 4px var(--gs)}to{box-shadow:0 0 16px var(--ge)}}
#alhdr span{color:var(--tc);font-size:11px;font-weight:bold;letter-spacing:2px;
  text-transform:uppercase;text-shadow:0 0 8px var(--gs);transition:color .35s}

/* ── Title ─────────────────────────────────────────────────────────────── */
h1{position:relative;z-index:1;color:var(--tc2);text-shadow:0 0 12px var(--gs);
  font-size:12.5px;letter-spacing:3.5px;text-transform:uppercase;text-align:center;
  margin-bottom:9px;line-height:1.4;transition:color .35s,text-shadow .35s}
h1 .tri{color:var(--tc);margin-right:6px}

/* ── Body ──────────────────────────────────────────────────────────────── */
#body{position:relative;z-index:1;display:flex;flex-direction:column;gap:6px}

/* ── Selects ───────────────────────────────────────────────────────────── */
.sr{display:flex;flex-direction:column;gap:2px}
.sr label{font-size:11px;letter-spacing:1.5px;text-transform:uppercase;color:var(--lc);transition:color .35s}
.sr select{background:var(--ibg);color:var(--tc2);border:1px solid var(--tc3);border-radius:2px;
  padding:5px 26px 5px 8px;font-family:'Courier New',monospace;font-size:12.5px;
  width:100%;cursor:pointer;outline:none;-webkit-appearance:none;appearance:none;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%2300e5cc'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 9px center;
  transition:background-color .35s,border-color .35s,color .35s}
.sr select:focus{border-color:var(--tc);box-shadow:0 0 6px var(--gs)}
.sr select option{background:#000d0d;color:#00e5cc}

/* ── Divider ───────────────────────────────────────────────────────────── */
.hdiv{border:none;border-top:1px solid var(--tc3);margin:2px 0;transition:border-color .35s}

/* ── Input fields ──────────────────────────────────────────────────────── */
#ifields{display:flex;flex-direction:column;gap:4px}
.frow{display:flex;flex-direction:column;gap:1px}
.frow label{font-size:11px;color:var(--lc);transition:color .35s}
.frow input{width:100%;background:var(--ibg);color:var(--tc2);
  border:1px solid var(--tc3);border-radius:2px;
  padding:5px 8px;font-family:'Courier New',monospace;font-size:12.5px;outline:none;
  transition:background-color .35s,border-color .35s,color .35s}
.frow input:focus{border-color:var(--tc);box-shadow:0 0 5px var(--gs)}
.frow input::-webkit-inner-spin-button,.frow input::-webkit-outer-spin-button{-webkit-appearance:none}
.frow input[type=number]{-moz-appearance:textfield}
.cap{font-size:10px;color:rgba(0,229,204,.32);letter-spacing:.5px}

/* ── Calculate button ──────────────────────────────────────────────────── */
#cbtn{width:100%;background:var(--bbg);color:var(--tc2);
  border:1.5px solid var(--tc);border-radius:2px;padding:7px 0;
  font-family:'Courier New',monospace;font-size:12.5px;letter-spacing:4px;
  cursor:pointer;text-transform:uppercase;
  transition:background .2s,box-shadow .2s,border-color .35s,color .35s}
#cbtn:hover{background:var(--bh);box-shadow:0 0 14px var(--ge)}
#cbtn:active{filter:brightness(1.3)}

/* ── Result boxes ──────────────────────────────────────────────────────── */
.rbox{padding:7px 11px;border-radius:2px;font-size:12.5px;line-height:1.65;word-break:break-word}
.ok{background:rgba(0,50,25,.93);border:1px solid rgba(0,190,95,.33);border-left:3px solid #00e87a;color:#00ffcc}
.ok .rt{font-size:13.5px;font-weight:bold}
.warn{background:rgba(44,26,0,.95);border:1px solid rgba(255,175,35,.30);border-left:3px solid #ffb800;color:#ffe590}
.blk{background:rgba(46,0,0,.97);border:1px solid rgba(255,55,55,.33);border-left:3px solid #ff3333;
  color:#ff8080;font-weight:bold;animation:bpls 1.5s ease-in-out infinite alternate}
@keyframes bpls{from{box-shadow:0 0 3px rgba(255,60,60,.16)}to{box-shadow:0 0 12px rgba(255,60,60,.46)}}
.err{background:rgba(35,0,0,.95);border:1px solid rgba(255,55,55,.28);border-left:3px solid #ff4444;color:#ff9090}

/* ── Spec card ─────────────────────────────────────────────────────────── */
.sc{margin-top:5px;padding:7px 11px;border-radius:2px;background:rgba(0,24,19,.97);
  border:1px solid var(--tc3);border-left:3px solid var(--tc);transition:border-color .35s}
.sh{color:rgba(0,229,204,.42);font-size:9px;letter-spacing:2px;text-transform:uppercase;margin-bottom:4px}
.sr2{color:#00ffee;font-size:12.5px;line-height:1.85}
.sn{color:rgba(0,229,204,.38);font-size:10px;margin-top:2px}
.badge{display:inline-block;border-radius:2px;font-size:8.5px;letter-spacing:1.5px;
  padding:1px 5px;margin-left:5px;vertical-align:middle;text-transform:uppercase}
.bo2{background:rgba(255,100,0,.13);border:1px solid rgba(255,120,0,.36);color:#ffaa44}
.bh2{background:rgba(0,175,255,.10);border:1px solid rgba(0,195,255,.30);color:#66ddff}

/* ── Confidence bar ────────────────────────────────────────────────────── */
#confbox{margin-top:5px;padding:5px 9px;border:1px solid var(--tc3);border-radius:2px;
  display:none;flex-direction:column;gap:3px;transition:border-color .35s}
#confbox.show{display:flex}
.clbl{font-size:8.5px;letter-spacing:2px;text-transform:uppercase;color:var(--lc)}
.ctrack{height:8px;background:rgba(0,0,0,.32);border-radius:1px;overflow:hidden;position:relative}
.cfill{height:100%;background:var(--tc);transition:width .6s ease,background .35s;position:relative}
.cfill::after{content:'';position:absolute;right:0;top:0;width:3px;height:100%;
  background:var(--tc2);opacity:.72;box-shadow:0 0 5px var(--gs)}
.cpct{font-size:10px;color:var(--tc);font-weight:bold;letter-spacing:1px;text-align:right;transition:color .35s}
</style>
</head>
<body>

%%HUD%%

<div id="wrapper">
  <div id="panel">

    <div id="alhdr"><span id="altxt">WARNING</span></div>

    <h1><span class="tri">&#9651;</span>HIGH-PRESSURE GAS FLOW CALCULATOR</h1>

    <div id="body">
      <div class="sr">
        <label>Gas Type</label>
        <select id="gasSelect" onchange="onGasChange()">
          <option value="N2">N2 (Nitrogen)</option>
          <option value="O2">O2 (Oxygen)</option>
          <option value="Ar">Ar (Argon)</option>
          <option value="CO2">CO2 (Carbon Dioxide)</option>
          <option value="He">He (Helium)</option>
          <option value="H2">H2 (Hydrogen)</option>
          <option value="CH4">CH4 (Methane)</option>
          <option value="C2H2">C2H2 (Acetylene)</option>
          <option value="FG1">(H2-3% + Ar-97%) Forming Gas 1</option>
          <option value="FG2">(H2-5% + N2-95%) Forming Gas 2</option>
          <option value="Air">Air (Dry Air)</option>
        </select>
      </div>

      <div class="sr">
        <label>Calculation Type</label>
        <select id="calcSelect" onchange="rebuildFields()">
          <option value="diameter">Pipe Diameter (mm)</option>
          <option value="flow">Flow Rate (LPM)</option>
          <option value="length">Pipe Length (m)</option>
          <option value="inlet">Inlet Pressure (bar)</option>
          <option value="outlet">Outlet Pressure (bar)</option>
        </select>
      </div>

      <hr class="hdiv">
      <div id="ifields"></div>
      <p class="cap">Friction factor (f) = 0.02 &nbsp;&middot;&nbsp; fixed constant</p>
      <button id="cbtn" onclick="doCalc()">&#9658;&nbsp;&nbsp;Calculate</button>
      <div id="rarea"></div>
      <div id="confbox">
        <div class="clbl">Recommendation Confidence</div>
        <div class="ctrack"><div class="cfill" id="cfill" style="width:0%"></div></div>
        <div class="cpct" id="cpct">0%</div>
      </div>
    </div>

  </div>
</div>

<script>
// TUBE_TABLE injected from Python
var TT=%%TUBE_TABLE%%;

// Safety block limits (inlet pressure → red panel + blocked calc)
var SL={
  "CO2" :{bar:69.9,title:"ALERT: SAFETY / LIMIT CONDITION",
          msg:"Risk: Liquid Phase CO\u2082\nPressure Exceeds Limit!",
          he:"\u05dc\u05d7\u05e5 \u05db\u05e0\u05d9\u05e1\u05d4 \u05d7\u05d5\u05e8\u05d2 \u2014 CO\u2082 \u05e2\u05d5\u05d1\u05e8 \u05dc\u05e0\u05d5\u05d6\u05dc"},
  "CH4" :{bar:250, title:"ALERT: SAFETY / LIMIT CONDITION",
          msg:"Risk: Flammable Gas\nPressure Exceeds Safe Limit!",
          he:"\u05dc\u05d7\u05e5 \u05db\u05e0\u05d9\u05e1\u05d4 \u05d7\u05d5\u05e8\u05d2 \u2014 \u05de\u05ea\u05d0\u05df \u05d1\u05dc\u05d7\u05e5 \u05de\u05e7\u05e1\u05d9\u05de\u05dc\u05d9"},
  "C2H2":{bar:17,  title:"ALERT: SAFETY / LIMIT CONDITION",
          msg:"Risk: Decomposition Hazard\nMax 17 bar for Acetylene!",
          he:"\u05d1\u05d8\u05d9\u05d7\u05d5\u05ea \u2014 \u05de\u05e7\u05e1 17 \u05d1\u05e8 \u05dc\u05d0\u05e6\u05d8\u05d9\u05dc\u05df"},
  "H2S" :{bar:20,  title:"ALERT: SAFETY / LIMIT CONDITION",
          msg:"Risk: Toxic Gas\nPressure Exceeds Safe Limit!",
          he:"\u05d1\u05d8\u05d9\u05d7\u05d5\u05ea \u2014 \u05dc\u05d7\u05e5 \u05de\u05e7\u05e1 20 \u05d1\u05e8 \u05dc-H\u2082S"}
};

// Gases that always show orange panel (dangerous regardless of pressure)
var WG=["O2","H2","CH4","C2H2","H2S"];

// Per-gas warning messages (orange state)
var WM={
  "O2" :"O\u2082 \u2014 Oxidizing gas. Use O\u2082-cleaned fittings. S14 spec only.",
  "H2" :"H\u2082 \u2014 Highly flammable. Ensure leak-free connections. S6 / S9 spec.",
  "CH4":"CH\u2084 \u2014 Flammable. Ensure adequate ventilation.",
  "C2H2":"C\u2082H\u2082 \u2014 Unstable above 17 bar. Keep below decomposition limit.",
  "H2S":"H\u2082S \u2014 Highly toxic. Gas detection and PPE mandatory."
};

// ── Theme engine ─────────────────────────────────────────────────────────────
function setTheme(m,title){
  document.getElementById("panel").className=m==="W"?"W":m==="D"?"D":"";
  var hdr=document.getElementById("alhdr"),txt=document.getElementById("altxt");
  hdr.style.display=(m==="N")?"none":"flex";
  if(txt&&title)txt.textContent=title;
}
function evalTheme(gas,P){
  var lim=SL[gas];
  if(lim&&P>lim.bar)return{m:"D",title:lim.title,blocked:true,lim:lim};
  if(WG.indexOf(gas)!==-1)return{m:"W",title:"WARNING: PRESSURE LIMIT REACHED",blocked:false};
  return{m:"N",blocked:false};
}

// ── CRITICAL FIX: toKey matches original — parens REMOVED (not replaced) ────
// "Inlet Pressure (bar)" -> "Inlet_Pressure_bar"
// "Temperature (°C)"     -> "Temperature_degC"
function toKey(f){
  return f.replace(/ /g,"_").replace(/[()]/g,"").replace(/\//g,"").replace(/\u00b0/g,"deg");
}

function getInletP(){
  var e=document.getElementById(toKey("Inlet Pressure (bar)"));
  return e?(parseFloat(e.value)||0):0;
}

function onGasChange(){
  rebuildFields();
  var gas=document.getElementById("gasSelect").value;
  var th=evalTheme(gas,getInletP());
  setTheme(th.m,th.title);
  if(th.m==="W"&&WM[gas]){
    document.getElementById("rarea").innerHTML=
      "<div class='rbox warn'><strong>\u26a0 "+gas+" \u2014 Hazardous Gas</strong><br>"+WM[gas]+"</div>";
  } else if(th.m==="N"){
    document.getElementById("rarea").innerHTML="";
  }
  hideConf();
}

function onPressureChange(){
  var gas=document.getElementById("gasSelect").value;
  setTheme(evalTheme(gas,getInletP()).m,evalTheme(gas,getInletP()).title);
}

// ── Tube lookup (unchanged from original) ────────────────────────────────────
function lookupTube(gas,inP,diam){
  var all=TT.filter(function(r){return r.gas_codes.indexOf(gas)!==-1&&r.max_pressure>=inP&&r.id_mm>=diam});
  if(!all.length)return null;
  var cap=all.filter(function(r){return r.max_pressure<=inP*5});
  var c=cap.length?cap:all;
  c.sort(function(a,b){return(a.id_mm-b.id_mm)||(a.max_pressure-b.max_pressure)});
  return c[0];
}
function fmtTube(s){return s.replace(/X/g," \u00d7 ").replace(/\s{2,}/g," ").trim()}
function buildSpec(gas,inP,diam){
  var m=lookupTube(gas,inP,diam);
  if(!m)return"<div class='rbox warn' style='margin-top:5px'><strong>\u26a0 No approved tube specification found</strong></div>";
  var badge="";
  if(gas==="O2")badge="<span class='badge bo2'>O\u2082 \u2192 S14 only</span>";
  if(gas==="H2")badge="<span class='badge bh2'>H\u2082 \u2192 S6/S9</span>";
  return"<div class='sc'><div class='sh'>\u25ba Recommended Tube Specification"+badge+"</div>"
    +"<div class='sr2'>\u2713&nbsp;<strong>Tube Spec:&nbsp;&nbsp;&nbsp;&nbsp;</strong>"+m.spec+"<br>"
    +"\u2713&nbsp;<strong>Min. Tube Size:&nbsp;</strong>"+fmtTube(m.tube_od_w)+"</div>"
    +"<div class='sn'>ID "+m.id_mm.toFixed(3)+" mm&nbsp;&nbsp;\u2502&nbsp;&nbsp;Max pressure rated: "+m.max_pressure+" bar</div></div>";
}

// Confidence score
function calcConf(gas,inP,diam,spec){
  if(!spec)return 5;
  var dr=Math.min((spec.id_mm-diam)/spec.id_mm,1);
  var pr=Math.max((spec.max_pressure-inP)/spec.max_pressure,0);
  var sf=WG.indexOf(gas)!==-1?0.72:1.0;
  return Math.max(5,Math.min(100,Math.round((dr*0.35+pr*0.65)*sf*100)));
}
function showConf(pct){
  var b=document.getElementById("confbox"),f=document.getElementById("cfill"),l=document.getElementById("cpct");
  b.className="show";setTimeout(function(){f.style.width=pct+"%";l.textContent=pct+"%"},60);
}
function hideConf(){var b=document.getElementById("confbox");if(b)b.className=""}

function renderResult(line,diam,inP,gas){
  var spec=lookupTube(gas,inP,diam);
  showConf(calcConf(gas,inP,diam,spec));
  return"<div class='rbox ok'><div class='rt'>"+line+"</div></div>"+buildSpec(gas,inP,diam);
}

// ── Gas physics (UNCHANGED — Math.pow instead of ** for compatibility) ────────
var GM={N2:.028013,O2:.031999,Ar:.039948,CO2:.04401,
        He:.0040026,H2:.002016,CH4:.01604,C2H2:.02604,
        FG1:.03881,FG2:.02671,Air:.02897};
var R=8.314,FR=0.02;
function rho(P,T,g){return(P*1e5*GM[g])/(R*(T+273.15))}
function rhoA(Pi,Po,T,g){return(rho(Pi,T,g)+rho(Po,T,g))/2}

function calcDiameter(Pi,Po,T,L,Q,g){
  var dP=(Pi-Po)*1e5;if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  var Qs=Q/60000;
  return Math.pow((FR*L*8*rhoA(Pi,Po,T,g)*Qs*Qs)/(Math.pow(Math.PI,2)*dP),0.2)*1000;
}
function calcFlow(Pi,Po,T,L,D,g){
  var dP=(Pi-Po)*1e5;if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  var Dm=D/1000;
  return Math.sqrt(dP*Math.pow(Math.PI,2)*Math.pow(Dm,5)/(8*FR*L*rhoA(Pi,Po,T,g)))*60000;
}
function calcLength(Pi,Po,T,D,Q,g){
  var dP=(Pi-Po)*1e5;if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  var Dm=D/1000,Qs=Q/60000;
  return dP*Math.pow(Math.PI,2)*Math.pow(Dm,5)/(8*FR*rhoA(Pi,Po,T,g)*Qs*Qs);
}
function calcOutlet(Pi,T,L,D,Q,g){
  var Dm=D/1000,Qs=Q/60000;
  function res(Po){return(Pi-Po)*1e5-(8*FR*L*rhoA(Pi,Po,T,g)*Qs*Qs)/(Math.pow(Math.PI,2)*Math.pow(Dm,5))}
  var lo=0,hi=Pi;
  for(var i=0;i<60;i++){if(Math.abs(hi-lo)<1e-4)break;var m=(lo+hi)/2;res(m)>0?lo=m:hi=m}
  return Math.max((lo+hi)/2,0);
}
function calcInlet(Po,T,L,D,Q,g){
  var lo=Po,hi=Po+10;
  while(hi<Po+2000){if(calcOutlet(hi,T,L,D,Q,g)>=Po)break;hi+=10}
  for(var i=0;i<60;i++){var m=(lo+hi)/2,vm=calcOutlet(m,T,L,D,Q,g);
    if(Math.abs(vm-Po)<0.005)return m;vm<Po?lo=m:hi=m}
  return(lo+hi)/2;
}

// ── Field configuration ──────────────────────────────────────────────────────
var FIELDS={
  diameter:["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  flow:    ["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  length:  ["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  inlet:   ["Temperature (\u00b0C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  outlet:  ["Temperature (\u00b0C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"]
};
var DEF={
  "Temperature (\u00b0C)":25,
  "Inlet Pressure (bar)":100,
  "Outlet Pressure (bar)":10,
  "Pipe Length (m)":10,
  "Pipe Diameter (mm)":10,
  "Flow Rate (LPM)":100
};
function fv(l){var e=document.getElementById(toKey(l));return e?(parseFloat(e.value)||0):(DEF[l]||0)}

// ── Main calculation ─────────────────────────────────────────────────────────
function doCalc(){
  var gas=document.getElementById("gasSelect").value;
  var ct =document.getElementById("calcSelect").value;
  var area=document.getElementById("rarea");
  try{
    var Tc=fv("Temperature (\u00b0C)"),Pi=fv("Inlet Pressure (bar)"),
        Po=fv("Outlet Pressure (bar)"),L=fv("Pipe Length (m)"),
        D=fv("Pipe Diameter (mm)"),Q=fv("Flow Rate (LPM)");

    var checkP=(ct==="inlet")?0:Pi;
    var th=evalTheme(gas,checkP);
    setTheme(th.m,th.title);

    // Red state: block calculation
    if(th.blocked&&ct!=="inlet"){
      var li=th.lim;
      area.innerHTML="<div class='rbox blk'>\u26a8 SAFETY BLOCKED<br>"
        +"<span style='font-weight:normal;font-size:11px'>"+li.msg.replace(/\n/g,"<br>")+"</span>"
        +"<div style='direction:rtl;font-size:10px;opacity:.76;margin-top:3px'>"+li.he+"</div></div>";
      showConf(5);return;
    }

    var html="";
    if(ct==="diameter"){
      var r=calcDiameter(Pi,Po,Tc,L,Q,gas);
      html=renderResult("Required Diameter: <strong>"+r.toFixed(2)+" mm</strong>",r,Pi,gas);
    }else if(ct==="flow"){
      var r=calcFlow(Pi,Po,Tc,L,D,gas);
      html=renderResult("Maximum Flow Rate: <strong>"+r.toFixed(1)+" L/min</strong>",D,Pi,gas);
    }else if(ct==="length"){
      var r=calcLength(Pi,Po,Tc,D,Q,gas);
      html=renderResult("Maximum Pipe Length: <strong>"+r.toFixed(1)+" m</strong>",D,Pi,gas);
    }else if(ct==="inlet"){
      var r=calcInlet(Po,Tc,L,D,Q,gas);
      var th2=evalTheme(gas,r);setTheme(th2.m,th2.title);
      if(th2.blocked){
        var li=th2.lim;
        html="<div class='rbox ok'><div class='rt'>Required Inlet Pressure: <strong>"+r.toFixed(2)+" bar</strong></div></div>"
          +"<div class='rbox blk' style='margin-top:5px'>\u26a8 Calculated pressure exceeds safe limit for "+gas+"<br>"
          +"<span style='font-weight:normal;font-size:11px'>"+li.msg.replace(/\n/g,"<br>")+"</span>"
          +"<div style='direction:rtl;font-size:10px;opacity:.76;margin-top:3px'>"+li.he+"</div></div>";
        showConf(10);
      }else{
        html=renderResult("Required Inlet Pressure: <strong>"+r.toFixed(2)+" bar</strong>",D,r,gas);
      }
    }else if(ct==="outlet"){
      var r=calcOutlet(Pi,Tc,L,D,Q,gas);
      html=renderResult("Estimated Outlet Pressure: <strong>"+r.toFixed(2)+" bar</strong>",D,Pi,gas);
    }
    area.innerHTML=html;
  }catch(e){area.innerHTML="<div class='rbox err'>\u26a0 "+e.message+"</div>";hideConf()}
}

// ── Field builder — uses corrected toKey ─────────────────────────────────────
function rebuildFields(){
  var ct=document.getElementById("calcSelect").value;
  var fl=FIELDS[ct]||[];
  var html="";
  fl.forEach(function(f){
    var k=toKey(f);
    var extra=(f==="Inlet Pressure (bar)")?" oninput='onPressureChange()'":"";
    html+="<div class='frow'><label>"+f+"</label>"
         +"<input type='number' id='"+k+"' value='"+(DEF[f]||0)+"' step='any'"+extra+"></div>";
  });
  document.getElementById("ifields").innerHTML=html;
  document.getElementById("rarea").innerHTML="";
  hideConf();
}

// Initialise
rebuildFields();
</script>
</body>
</html>"""

page_html = PAGE.replace("%%HUD%%", HUD_TAG).replace("%%TUBE_TABLE%%", TUBE_TABLE_JSON)
components.html(page_html, height=950, scrolling=False)
