"""
High-Pressure Gas Flow Calculator
──────────────────────────────────
• HUD background (hud_background.html) fills the full screen and keeps animating
• Calculator panel sits EXACTLY inside the large black rectangle in the HUD
  (left:1%, top:21%, width:54%, height:44% – matches the empty rectangle visible
   in the screenshot)
• Panel has internal scroll (overflow-y:auto) – it NEVER overflows onto the HUD
• All physics runs in JavaScript – no page reload needed
• Border colour: #00e5cc (cyan/teal matching the HUD palette)

Run:  streamlit run app.py
Req:  hud_background.html in the same folder
"""

import base64
import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path

st.set_page_config(
    page_title="High-Pressure Gas Flow Calculator",
    page_icon="💨",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# Hide all Streamlit chrome; expand the components iframe to cover 100 % of viewport
st.markdown("""
<style>
  #MainMenu,footer,header,
  [data-testid="stToolbar"],[data-testid="stDecoration"],
  [data-testid="stStatusWidget"],[data-testid="collapsedControl"]{display:none!important}
  .stApp{background:#000d0d!important}
  section.main .block-container{padding:0!important;margin:0!important;max-width:100vw!important}
  iframe[title="streamlit_components_v1.html_v1"]{
    position:fixed!important;inset:0!important;
    width:100vw!important;height:100vh!important;
    border:none!important;z-index:100!important}
</style>
""", unsafe_allow_html=True)


def load_hud() -> str:
    p = Path(__file__).resolve().parent / "hud_background.html"
    if not p.exists():
        return ""
    return "data:text/html;base64," + base64.b64encode(p.read_bytes()).decode()


hud_src  = load_hud()
hud_tag  = (
    f"<iframe id='hud' src='{hud_src}' sandbox='allow-scripts'></iframe>"
    if hud_src else
    "<div style='position:fixed;inset:0;background:#000d0d;z-index:1'></div>"
)

# ─────────────────────────────────────────────────────────────────────────────
#  Full self-contained HTML page
#  The calculator panel is positioned to sit INSIDE the large empty rectangle
#  that occupies roughly left:1%, top:21%, width:54%, height:44% of the HUD.
# ─────────────────────────────────────────────────────────────────────────────
HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
html,body{
  width:100%;height:100%;
  overflow:hidden;          /* page itself NEVER scrolls */
  background:#000d0d;
  font-family:'Courier New',monospace;
  color:#00e5cc;
}

/* ── HUD: true full-screen, behind everything ── */
#hud{
  position:fixed;inset:0;
  width:100%;height:100%;
  border:none;pointer-events:none;
  z-index:1;
}

/* ══════════════════════════════════════════════════════════
   CALCULATOR PANEL
   Positioned to sit exactly inside the large black rectangle
   that is visible in the HUD screenshot:
     left ≈ 1%,  top ≈ 21%,  width ≈ 54%,  height ≈ 44%
   ══════════════════════════════════════════════════════════ */
#calc-panel{
  position:fixed;
  left:1%;
  top:21%;
  width:54%;
  height:44%;

  /* cyan border matching HUD palette */
  border:1.5px solid #00e5cc;
  border-radius:3px;
  box-shadow:0 0 18px rgba(0,229,204,0.35), inset 0 0 40px rgba(0,229,204,0.04);
  animation:panelGlow 2.5s ease-in-out infinite alternate;

  background:rgba(0,8,8,0.94);
  z-index:500;

  /* internal scroll ONLY */
  overflow-y:auto;
  overflow-x:hidden;
  padding:14px 18px 16px;

  scrollbar-width:thin;
  scrollbar-color:#00e5cc rgba(0,229,204,0.08);
}
#calc-panel::-webkit-scrollbar{width:5px}
#calc-panel::-webkit-scrollbar-track{background:rgba(0,229,204,0.05)}
#calc-panel::-webkit-scrollbar-thumb{background:#00e5cc;border-radius:3px}

@keyframes panelGlow{
  from{box-shadow:0 0 8px rgba(0,229,204,0.25),inset 0 0 20px rgba(0,229,204,0.03)}
  to  {box-shadow:0 0 24px rgba(0,229,204,0.6),inset 0 0 50px rgba(0,229,204,0.07)}
}

/* scanline overlay inside panel */
#calc-panel::after{
  content:"";position:absolute;inset:0;pointer-events:none;
  background:repeating-linear-gradient(
    to bottom,
    transparent 0px,transparent 3px,
    rgba(0,229,204,0.01) 3px,rgba(0,229,204,0.01) 4px);
  border-radius:3px;z-index:0;
}

/* ── Title ── */
#calc-panel h1{
  color:#00ffee;
  text-shadow:0 0 12px #00e5cc;
  font-size:clamp(11px,1.1vw,16px);
  letter-spacing:3px;text-transform:uppercase;
  text-align:center;margin-bottom:12px;line-height:1.4;
}
#calc-panel h1 .tri{color:#00e5cc;margin-right:6px}

/* ── Two-column layout inside panel ── */
#panel-body{
  display:grid;
  grid-template-columns:1fr 1fr;
  gap:10px 16px;
  position:relative;z-index:1;
}

/* ── Select rows ── */
.sr{display:flex;flex-direction:column;gap:3px;margin-bottom:8px}
.sr label{font-size:9px;letter-spacing:1.5px;text-transform:uppercase;color:rgba(0,229,204,0.6)}
.sr select{
  background:rgba(0,15,13,0.98);color:#00ffee;
  border:1px solid rgba(0,229,204,0.4);border-radius:2px;
  padding:5px 22px 5px 7px;
  font-family:'Courier New',monospace;font-size:10.5px;
  width:100%;cursor:pointer;outline:none;
  -webkit-appearance:none;appearance:none;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='9' height='5'%3E%3Cpath d='M0 0l4.5 5 4.5-5z' fill='%2300e5cc'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 8px center;
}
.sr select:focus{border-color:#00e5cc;box-shadow:0 0 6px rgba(0,229,204,0.35)}
.sr select option{background:#000d0d;color:#00e5cc}

/* top selects span both columns */
.span2{grid-column:1/-1}

.hdiv{grid-column:1/-1;border:none;border-top:1px solid rgba(0,229,204,0.18);margin:2px 0 4px}

/* ── Input field rows ── */
.frow{display:flex;flex-direction:column;gap:3px}
.frow label{font-size:9px;color:rgba(0,229,204,0.75);letter-spacing:0.3px}
.frow input{
  width:100%;
  background:rgba(0,15,13,0.98);color:#00ffee;
  border:1px solid rgba(0,229,204,0.38);border-radius:2px;
  padding:5px 7px;font-family:'Courier New',monospace;font-size:11px;outline:none;
}
.frow input:focus{border-color:#00e5cc;box-shadow:0 0 5px rgba(0,229,204,0.35)}
.frow input::-webkit-inner-spin-button,
.frow input::-webkit-outer-spin-button{-webkit-appearance:none;margin:0}
.frow input[type=number]{-moz-appearance:textfield}

/* ── Caption + button ── */
.cap{grid-column:1/-1;font-size:8px;color:rgba(0,229,204,0.38);letter-spacing:1px}
#cbtn{
  grid-column:1/-1;
  background:rgba(0,229,204,0.07);color:#00ffee;
  border:1.5px solid #00e5cc;border-radius:2px;
  padding:7px 0;font-family:'Courier New',monospace;
  font-size:11px;letter-spacing:3px;cursor:pointer;text-transform:uppercase;
  transition:background .2s,box-shadow .2s;
}
#cbtn:hover{background:rgba(0,229,204,0.18);box-shadow:0 0 14px rgba(0,229,204,0.5)}
#cbtn:active{background:rgba(0,229,204,0.28)}

/* ── Result box ── */
#result-area{grid-column:1/-1}
.rbox{
  padding:9px 12px;border-radius:2px;
  font-size:11px;line-height:1.85;word-break:break-word;
  margin-top:6px;
}
.rbox.ok{
  background:rgba(0,60,32,0.92);
  border:1px solid rgba(0,200,100,0.4);
  border-left:3px solid #00e87a;
  color:#00ffcc;font-weight:bold;
}
.rbox.ok small{font-weight:normal;font-size:10px;color:rgba(0,255,200,0.7)}
.rbox.err{
  background:rgba(40,0,0,0.94);
  border:1px solid rgba(255,60,60,0.3);
  border-left:3px solid #ff4444;
  color:#ff9090;
}

/* ── Spec table ── */
.sp-title{font-size:9px;letter-spacing:1px;color:#00e5cc;text-transform:uppercase;margin:8px 0 4px;font-weight:bold}
.stbl{width:100%;border-collapse:collapse;font-size:10px}
.stbl th{
  color:#00e5cc;background:rgba(0,229,204,0.07);
  border-bottom:1px solid rgba(0,229,204,0.28);
  padding:4px 7px;text-align:left;letter-spacing:1px;font-size:9px;text-transform:uppercase;
}
.stbl td{
  color:rgba(0,229,204,0.82);border-bottom:1px solid rgba(0,229,204,0.1);padding:4px 7px;
}
.stbl tr:last-child td{border-bottom:none}
</style>
</head>
<body>

%%HUD%%

<div id="calc-panel">
  <h1><span class="tri">&#9651;</span>HIGH-PRESSURE GAS FLOW CALCULATOR</h1>

  <div id="panel-body">

    <div class="sr span2">
      <label>Gas type</label>
      <select id="gasSelect" onchange="rebuildFields()">
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

    <div class="sr span2">
      <label>Calculation type</label>
      <select id="calcSelect" onchange="rebuildFields()">
        <option value="diameter">Pipe Diameter (mm)</option>
        <option value="flow">Flow Rate (LPM)</option>
        <option value="length">Pipe Length (m)</option>
        <option value="inlet">Inlet Pressure (bar)</option>
        <option value="outlet">Outlet Pressure (bar)</option>
      </select>
    </div>

    <hr class="hdiv">

    <div id="input-fields"></div>

    <p class="cap">Friction factor (f) = 0.02 &nbsp;&middot;&nbsp; fixed constant</p>
    <button id="cbtn" onclick="doCalc()">&#9658;&nbsp;&nbsp;Calculate</button>

    <div id="result-area"></div>

  </div><!-- #panel-body -->
</div><!-- #calc-panel -->

<script>
// ── Gas molecular weights kg/mol ──────────────────────────────────────────
const GAS_M={N2:.028013,O2:.031999,Ar:.039948,CO2:.04401,He:.0040026,
             H2:.002016,CH4:.01604,C2H2:.02604,FG1:.03881,FG2:.02671,Air:.02897};
const R=8.314, FR=0.02;

// ── Field definitions per calc type ──────────────────────────────────────
const FIELDS={
  diameter:["\u0054emperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  flow:    ["\u0054emperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  length:  ["\u0054emperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  inlet:   ["\u0054emperature (\u00b0C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  outlet:  ["\u0054emperature (\u00b0C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"]
};
const DEF={
  "\u0054emperature (\u00b0C)":25,"Inlet Pressure (bar)":100,
  "Outlet Pressure (bar)":10,"Pipe Length (m)":10,
  "Pipe Diameter (mm)":10,"Flow Rate (LPM)":100
};

function toKey(f){return f.replace(/ /g,"_").replace(/[()]/g,"").replace(/\//g,"").replace(/\u00b0/g,"deg").replace(/\u0054/g,"T")}
function fv(label){const el=document.getElementById(toKey(label));return el?parseFloat(el.value)||0:(DEF[label]||0)}

// ── Physics ───────────────────────────────────────────────────────────────
function rho(Pbar,Tc,gas){return(Pbar*1e5*GAS_M[gas])/(R*(Tc+273.15))}
function rhoAvg(Pi,Po,Tc,gas){return(rho(Pi,Tc,gas)+rho(Po,Tc,gas))/2}

function calcDiameter(Pi,Po,Tc,L,Q,gas){
  const dP=(Pi-Po)*1e5;
  if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  return Math.pow((FR*L*8*rhoAvg(Pi,Po,Tc,gas)*(Q/60000)**2)/(Math.PI**2*dP),0.2)*1000;
}
function calcFlow(Pi,Po,Tc,L,D,gas){
  const dP=(Pi-Po)*1e5;
  if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  return Math.sqrt(dP*Math.PI**2*(D/1000)**5/(8*FR*L*rhoAvg(Pi,Po,Tc,gas)))*60000;
}
function calcLength(Pi,Po,Tc,D,Q,gas){
  const dP=(Pi-Po)*1e5;
  if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  return dP*Math.PI**2*(D/1000)**5/(8*FR*rhoAvg(Pi,Po,Tc,gas)*(Q/60000)**2);
}
function calcOutlet(Pi,Tc,L,D,Q,gas){
  const Dm=D/1000,Qs=Q/60000;
  function res(Po){return(Pi-Po)*1e5-(8*FR*L*rhoAvg(Pi,Po,Tc,gas)*Qs**2)/(Math.PI**2*Dm**5)}
  let lo=0,hi=Pi;
  for(let i=0;i<60;i++){if(Math.abs(hi-lo)<1e-4)break;const m=(lo+hi)/2;res(m)>0?lo=m:hi=m}
  return Math.max((lo+hi)/2,0);
}
function calcInlet(Po,Tc,L,D,Q,gas){
  let lo=Po,hi=Po+10;
  while(hi<Po+2000){if(calcOutlet(hi,Tc,L,D,Q,gas)>=Po)break;hi+=10}
  for(let i=0;i<60;i++){
    const m=(lo+hi)/2,vm=calcOutlet(m,Tc,L,D,Q,gas);
    if(Math.abs(vm-Po)<.005)return m;
    vm<Po?lo=m:hi=m;
  }
  return(lo+hi)/2;
}

// ── Pipe spec lookup ──────────────────────────────────────────────────────
function pipeSpec(D,P,gas){
  if(gas==="O2")return{rec:'1" (S14)',rows:[['1" (S14)','Required for O2']]};
  let rec="",rows=[];
  if(D<=4){
    if(P<=200){rows=[['1/4" (S6)','&#8804;200 bar'],['1/4" (S9)','&#8804;1379 bar'],['1/4" (S12)','>1379 bar']];rec='1/4" (S6)'}
    else if(P<=1379){rows=[['1/4" (S9)','&#8804;1379 bar'],['1/4" (S12)','>1379 bar']];rec='1/4" (S9)'}
    else{rows=[['1/4" (S12)','>1379 bar']];rec='1/4" (S12)'}
  }else if(D<=7){
    if(P<=140){rows=[['3/8" (S16)','&#8804;140 bar'],['3/8" (S9)','>140 bar']];rec='3/8" (S16)'}
    else{rows=[['3/8" (S9)','>140 bar']];rec='3/8" (S9)'}
  }else if(D<=21){
    if(P<=20){rows=[['3/4" (S15)','&#8804;20 bar'],['1" (S14)','>20 bar']];rec='3/4" (S15)'}
    else{rows=[['1" (S14)','>20 bar']];rec='1" (S14)'}
  }else{rows=[['Special piping','Outside standard range']];rec='Special'}
  return{rec,rows};
}

function renderResult(mainLine,D,P,gas){
  const{rec,rows}=pipeSpec(D,P,gas);
  const tr=rows.map(r=>"<tr><td>"+r[0]+"</td><td>"+r[1]+"</td></tr>").join("");
  return "<div class='rbox ok'>"+mainLine+"<br><small>Recommended: "+rec+"</small></div>"
        +"<p class='sp-title'>Possible pipe specifications:</p>"
        +"<table class='stbl'><thead><tr><th>Pipe Spec</th><th>Range</th></tr></thead>"
        +"<tbody>"+tr+"</tbody></table>";
}

// ── CALCULATE ─────────────────────────────────────────────────────────────
function doCalc(){
  const gas=document.getElementById("gasSelect").value;
  const ct =document.getElementById("calcSelect").value;
  const area=document.getElementById("result-area");
  try{
    const Tc=fv("Temperature (\u00b0C)"),Pi=fv("Inlet Pressure (bar)"),
          Po=fv("Outlet Pressure (bar)"),L=fv("Pipe Length (m)"),
          D=fv("Pipe Diameter (mm)"),Q=fv("Flow Rate (LPM)");
    let html="";
    if(ct==="diameter"){const r=calcDiameter(Pi,Po,Tc,L,Q,gas);html=renderResult("Required Diameter: <strong>"+r.toFixed(2)+" mm</strong>",r,Math.max(Pi,Po),gas)}
    else if(ct==="flow"){const r=calcFlow(Pi,Po,Tc,L,D,gas);html=renderResult("Maximum Flow Rate: <strong>"+r.toFixed(1)+" L/min</strong>",D,Math.max(Pi,Po),gas)}
    else if(ct==="length"){const r=calcLength(Pi,Po,Tc,D,Q,gas);html=renderResult("Maximum Pipe Length: <strong>"+r.toFixed(1)+" m</strong>",D,Math.max(Pi,Po),gas)}
    else if(ct==="inlet"){const r=calcInlet(Po,Tc,L,D,Q,gas);html=renderResult("Required Inlet Pressure: <strong>"+r.toFixed(2)+" bar</strong>",D,Math.max(r,Po),gas)}
    else if(ct==="outlet"){const r=calcOutlet(Pi,Tc,L,D,Q,gas);html=renderResult("Estimated Outlet Pressure: <strong>"+r.toFixed(2)+" bar</strong>",D,Pi,gas)}
    area.innerHTML=html;
    // scroll result into view
    area.scrollIntoView({behavior:"smooth",block:"nearest"});
  }catch(e){area.innerHTML="<div class='rbox err'>&#9888; "+e.message+"</div>"}
}

// ── Rebuild fields when calc type changes ─────────────────────────────────
function toKeyLocal(f){return f.replace(/ /g,"_").replace(/[()]/g,"").replace(/\//g,"").replace(/\u00b0/g,"deg")}

function rebuildFields(){
  const ct=document.getElementById("calcSelect").value;
  const fl=FIELDS[ct]||[];
  // split into two columns (pairs)
  let html="";
  fl.forEach(function(f){
    const k=toKeyLocal(f);
    html+="<div class='frow'><label>"+f+"</label>"
         +"<input type='number' id='"+k+"' value='"+(DEF[f]||0)+"' step='any'></div>";
  });
  document.getElementById("input-fields").innerHTML=html;
  document.getElementById("result-area").innerHTML="";
}

// Make input-fields span 2 cols inside the grid
document.getElementById("input-fields").style.gridColumn="1 / -1";
document.getElementById("result-area").style.gridColumn="1 / -1";

rebuildFields();
</script>
</body>
</html>"""

page_html = HTML.replace("%%HUD%%", hud_tag)
components.html(page_html, height=800, scrolling=False)
