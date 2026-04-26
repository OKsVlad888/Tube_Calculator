"""
High-Pressure Gas Flow Calculator  —  v4
=========================================
Key changes vs v3:
  - Panel positioned to sit EXACTLY inside the HUD black center zone
    (left 22%, right 39.7%, top 5.6% — matches HUD cleared area)
  - Panel fills the full black zone (no fixed max-width)
  - overflow-y:auto so it scrolls inside its zone only
  - All theme tokens (teal / orange / red) apply to every element
  - Confidence bar shown in ALL modes (normal included)
  - 2-column field layout in warn/danger modes only

Run: streamlit run app.py
Req: hud_background.html in the same directory
"""

import base64, json
import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path

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
  .stApp{background:#010a08!important}
  section.main,.block-container{padding:0!important;margin:0!important;max-width:100vw!important}
  iframe[title="streamlit_components_v1.html_v1"]{
    position:fixed!important;inset:0!important;
    width:100vw!important;height:100vh!important;
    border:none!important;z-index:100!important;
    display:block!important}
</style>
""", unsafe_allow_html=True)

# ── TUBE TABLE (unchanged) ────────────────────────────────────────────────────
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

def load_hud_tag() -> str:
    p = Path(__file__).resolve().parent / "hud_background.html"
    if p.exists():
        b64 = base64.b64encode(p.read_bytes()).decode()
        return (
            "<iframe src='data:text/html;base64," + b64 + "'"
            " sandbox='allow-scripts'"
            " style='position:fixed;inset:0;width:100%;height:100%;"
            "border:none;z-index:1;pointer-events:none'></iframe>"
        )
    return "<div style='position:fixed;inset:0;background:#010a08;z-index:1'></div>"

HUD_TAG = load_hud_tag()

PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
html,body{width:100%;height:100%;overflow:hidden;background:#010a08;
  font-family:'Courier New',monospace}

/* ── CSS theme variables ─────────────────────────────────────────────────── */
:root{
  --tc:#00e5cc;  --tc2:#00ffee; --tc3:rgba(0,229,204,0.38);
  --tc4:rgba(0,229,204,0.20);  --lc:rgba(0,229,204,0.78);
  --pbg:rgba(0,6,6,0.97);      --ibg:rgba(0,12,10,0.99);
  --bbg:rgba(0,229,204,0.07);  --bh:rgba(0,229,204,0.22);
  --gs:rgba(0,229,204,0.30);   --ge:rgba(0,229,204,0.75);
  --rbg:rgba(0,229,204,0.05);  --resbox-bg:rgba(0,28,22,0.98);
}
.warn{
  --tc:#ff8c00;  --tc2:#ffcc66; --tc3:rgba(255,140,0,0.50);
  --tc4:rgba(255,140,0,0.25);  --lc:rgba(255,190,60,0.88);
  --pbg:rgba(18,7,0,0.98);     --ibg:rgba(26,10,0,0.99);
  --bbg:rgba(255,140,0,0.07);  --bh:rgba(255,140,0,0.24);
  --gs:rgba(255,140,0,0.32);   --ge:rgba(255,140,0,0.82);
  --rbg:rgba(255,140,0,0.05);  --resbox-bg:rgba(26,10,0,0.98);
}
.dngr{
  --tc:#ff2222;  --tc2:#ff7777; --tc3:rgba(255,40,40,0.50);
  --tc4:rgba(255,40,40,0.25);  --lc:rgba(255,120,120,0.88);
  --pbg:rgba(18,0,0,0.99);     --ibg:rgba(24,2,2,0.99);
  --bbg:rgba(255,40,40,0.07);  --bh:rgba(255,40,40,0.24);
  --gs:rgba(255,40,40,0.32);   --ge:rgba(255,40,40,0.88);
  --rbg:rgba(255,40,40,0.05);  --resbox-bg:rgba(24,2,2,0.98);
}

/* ── Panel wrapper — sits exactly in HUD black zone ─────────────────────── */
/* HUD clears px(422)/1920 → px(1158)/1920, top py(58)/1037            */
/* CSS equivalents: left=22%, right=39.7%, top=5.6%                    */
#wrapper{
  position:fixed;
  left:22%;       /* matches HUD black zone left  (px422/1920) */
  right:39.7%;    /* matches HUD black zone right (1-px1158/1920) */
  top:5.6%;       /* matches HUD top bar          (py58/1037) */
  bottom:0;
  display:flex;flex-direction:column;
  z-index:500;pointer-events:none;
  padding:8px 10px 6px;
}

/* ── Panel fills the zone ────────────────────────────────────────────────── */
#panel{
  pointer-events:all;
  width:100%;flex:1;
  overflow-y:auto;overflow-x:hidden;
  background:var(--pbg);
  border:1.5px solid var(--tc);border-radius:4px;
  padding:12px 20px 14px;position:relative;
  animation:pglow 2.5s ease-in-out infinite alternate;
  transition:background 0.35s,border-color 0.35s;
  scrollbar-width:thin;scrollbar-color:var(--tc) transparent;
}
#panel::-webkit-scrollbar{width:4px}
#panel::-webkit-scrollbar-thumb{background:var(--tc);border-radius:2px}
@keyframes pglow{
  from{box-shadow:0 0 10px var(--gs),inset 0 0 28px rgba(0,0,0,0.5)}
  to  {box-shadow:0 0 38px var(--ge),inset 0 0 55px rgba(0,0,0,0.25)}}
#panel::after{content:"";position:absolute;inset:0;pointer-events:none;
  border-radius:4px;z-index:0;
  background:repeating-linear-gradient(to bottom,transparent 0,transparent 3px,
    rgba(0,229,204,0.007) 3px,rgba(0,229,204,0.007) 4px)}

/* ── Alert header ────────────────────────────────────────────────────────── */
#alert-hdr{
  display:none;align-items:center;justify-content:center;gap:9px;
  padding:7px 12px;margin-bottom:9px;
  background:var(--rbg);border:1px solid var(--tc3);border-radius:3px;
  animation:ahpls 1.2s ease-in-out infinite alternate;
  position:relative;z-index:1;transition:border-color 0.35s}
.warn #alert-hdr,.dngr #alert-hdr{display:flex}
@keyframes ahpls{from{box-shadow:0 0 4px var(--gs)}to{box-shadow:0 0 22px var(--ge)}}
.ah-icon{font-size:16px;color:var(--tc)}
.ah-text{font-size:13px;font-weight:bold;letter-spacing:2.5px;
  text-transform:uppercase;color:var(--tc);text-shadow:0 0 10px var(--gs)}

/* ── Title ───────────────────────────────────────────────────────────────── */
h1{position:relative;z-index:1;color:var(--tc2);text-shadow:0 0 14px var(--gs);
  font-size:15px;letter-spacing:4px;text-transform:uppercase;
  text-align:center;margin-bottom:10px;line-height:1.4;
  transition:color 0.35s,text-shadow 0.35s}
.tri{color:var(--tc);margin-right:7px}

/* ── Body ────────────────────────────────────────────────────────────────── */
#body{position:relative;z-index:1;display:flex;flex-direction:column;gap:8px}

/* ── Selects — normal stacked, 2-col in warn/dngr ───────────────────────── */
.sr{display:flex;flex-direction:column;gap:3px}
.sr label{font-size:11px;letter-spacing:1.5px;text-transform:uppercase;
  color:var(--lc);transition:color 0.35s}
.sr select{
  background:var(--ibg);color:var(--tc2);border:1px solid var(--tc3);
  border-radius:2px;padding:5px 28px 5px 9px;
  font-family:'Courier New',monospace;font-size:14px;
  width:100%;cursor:pointer;outline:none;
  -webkit-appearance:none;appearance:none;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%2300e5cc'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 10px center;
  transition:background-color 0.35s,border-color 0.35s,color 0.35s}
.sr select:focus{border-color:var(--tc);box-shadow:0 0 7px var(--gs)}
.sr select option{background:#000d0d;color:#00e5cc}
/* 2-col layout in non-normal mode */
.warn .sr,.dngr .sr{flex-direction:row;align-items:center;gap:10px}
.warn .sr label,.dngr .sr label{flex:0 0 44%;font-size:12px}

/* ── Divider ─────────────────────────────────────────────────────────────── */
.hdiv{border:none;border-top:1px solid var(--tc4);margin:1px 0;transition:border-color 0.35s}

/* ── Input fields ────────────────────────────────────────────────────────── */
#ifields{display:flex;flex-direction:column;gap:5px}
.frow{display:flex;flex-direction:column;gap:2px}
.frow label{font-size:13px;color:var(--lc);transition:color 0.35s}
.frow input{
  width:100%;background:var(--ibg);color:var(--tc2);
  border:1px solid var(--tc3);border-radius:2px;
  padding:5px 9px;font-family:'Courier New',monospace;font-size:14px;outline:none;
  transition:background-color 0.35s,border-color 0.35s,color 0.35s}
.frow input:focus{border-color:var(--tc);box-shadow:0 0 6px var(--gs)}
.frow input::-webkit-inner-spin-button,.frow input::-webkit-outer-spin-button{-webkit-appearance:none}
.frow input[type=number]{-moz-appearance:textfield}
/* 2-col in warn/dngr */
.warn .frow,.dngr .frow{flex-direction:row;align-items:center;gap:10px}
.warn .frow label,.dngr .frow label{flex:0 0 44%}
.warn .frow input,.dngr .frow input{flex:1}

/* ── Caption ─────────────────────────────────────────────────────────────── */
.cap{font-size:11px;color:var(--lc);opacity:0.45;letter-spacing:0.8px}

/* ── Calculate button ────────────────────────────────────────────────────── */
#cbtn{
  width:100%;background:var(--bbg);color:var(--tc2);
  border:1.5px solid var(--tc);border-radius:2px;padding:8px 0;
  font-family:'Courier New',monospace;font-size:14px;
  letter-spacing:4px;cursor:pointer;text-transform:uppercase;
  transition:background 0.2s,box-shadow 0.2s,border-color 0.35s,color 0.35s}
#cbtn:hover{background:var(--bh);box-shadow:0 0 20px var(--ge)}
#cbtn:active{filter:brightness(1.3)}
.dngr #cbtn{opacity:0.4;cursor:not-allowed;pointer-events:none}

/* ── Result boxes ────────────────────────────────────────────────────────── */
.rbox{padding:9px 13px;border-radius:2px;font-size:13px;line-height:1.7;word-break:break-word}
.rbox.ok{background:rgba(0,50,26,0.96);border:1px solid rgba(0,200,100,0.42);
  border-left:3px solid #00e87a;color:#00ffcc}
.rbox.ok .res-title{font-size:15px;font-weight:bold}
.rbox.blk{background:rgba(46,0,0,0.98);border:1px solid rgba(255,40,40,0.42);
  border-left:3px solid #ff2222;color:#ff7777;font-weight:bold;
  animation:bpls 1.4s ease-in-out infinite alternate}
@keyframes bpls{from{box-shadow:0 0 4px rgba(255,40,40,0.2)}to{box-shadow:0 0 16px rgba(255,40,40,0.55)}}
.rbox.err{background:rgba(34,0,0,0.96);border:1px solid rgba(255,50,50,0.30);
  border-left:3px solid #ff4444;color:#ff9090}

/* ── Spec card (normal mode) ─────────────────────────────────────────────── */
.spec-card{margin-top:7px;padding:9px 13px;border-radius:2px;
  background:rgba(0,24,20,0.98);border:1px solid var(--tc4);
  border-left:3px solid var(--tc);transition:border-color 0.35s}
.sc-head{color:var(--lc);opacity:0.55;font-size:10px;letter-spacing:2px;
  text-transform:uppercase;margin-bottom:5px}
.sc-row{color:var(--tc2);font-size:14px;line-height:1.95}
.sc-note{color:var(--lc);opacity:0.48;font-size:11px;margin-top:2px}
.badge{display:inline-block;border-radius:2px;font-size:9px;letter-spacing:1.5px;
  padding:1px 6px;margin-left:6px;vertical-align:middle;text-transform:uppercase}
.b-o2{background:rgba(255,100,0,0.15);border:1px solid rgba(255,120,0,0.4);color:#ffaa44}
.b-h2{background:rgba(0,180,255,0.12);border:1px solid rgba(0,200,255,0.35);color:#66ddff}

/* ── RESULTS box (warn/dngr mode) ───────────────────────────────────────── */
#res-box{
  display:none;padding:10px 13px;border-radius:2px;
  background:var(--resbox-bg);border:1px solid var(--tc4);
  transition:background 0.35s,border-color 0.35s}
.warn #res-box,.dngr #res-box{display:block}
.res-lbl{font-size:10px;letter-spacing:2.5px;text-transform:uppercase;
  color:var(--lc);opacity:0.75;margin-bottom:6px}
#res-content{font-size:13px;color:var(--tc2);line-height:1.85;transition:color 0.35s}

/* ── Confidence bar ──────────────────────────────────────────────────────── */
#conf-wrap{margin-top:5px;display:none}
#conf-wrap.show{display:block}
.conf-row{display:flex;align-items:center;gap:10px}
.conf-lbl{font-size:9px;letter-spacing:2px;text-transform:uppercase;
  color:var(--lc);white-space:nowrap;opacity:0.8}
.conf-track{
  flex:1;height:14px;background:rgba(0,0,0,0.38);
  border:1px solid var(--tc4);border-radius:1px;
  overflow:hidden;position:relative}
.conf-fill{
  height:100%;background:var(--tc);
  position:absolute;left:0;top:0;
  transition:width 0.75s ease,background 0.35s;
  display:flex;align-items:center;justify-content:flex-end;padding-right:3px}
/* Segmented notches */
.conf-fill::before{content:'';position:absolute;inset:0;
  background:repeating-linear-gradient(90deg,
    transparent 0,transparent calc(10% - 1px),
    rgba(0,0,0,0.40) calc(10% - 1px),rgba(0,0,0,0.40) 10%)}
.conf-pct{font-size:13px;font-weight:bold;color:var(--tc);
  min-width:42px;text-align:right;transition:color 0.35s}
</style>
</head>
<body>

%%HUD%%

<div id="wrapper">
<div id="panel">

  <div id="alert-hdr">
    <span class="ah-icon">&#9888;</span>
    <span class="ah-text" id="alert-txt">WARNING</span>
  </div>

  <h1><span class="tri">&#9651;</span>HIGH-PRESSURE GAS FLOW CALCULATOR</h1>

  <div id="body">

    <div class="sr">
      <label>Gas Type</label>
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

    <!-- Normal mode result -->
    <div id="result-normal"></div>

    <!-- Warn/Danger RESULTS box -->
    <div id="res-box">
      <div class="res-lbl">Results</div>
      <div id="res-content"></div>
    </div>

    <!-- Confidence bar (all modes) -->
    <div id="conf-wrap">
      <div class="conf-row">
        <span class="conf-lbl">Recommendation Confidence</span>
        <div class="conf-track"><div class="conf-fill" id="conf-fill" style="width:0%"></div></div>
        <span class="conf-pct" id="conf-pct">0%</span>
      </div>
    </div>

  </div>
</div>
</div>

<script>
var TT = %%TUBE_TABLE%%;

// ── Safety rules ──────────────────────────────────────────────────────────────
var BLK = {
  "CO2" :{max:69.9, hdr:"ALERT: SAFETY / LIMIT CONDITION",
    msg:"Risk: Liquid Phase CO\u2082\nPressure Exceeds Limit!"},
  "CH4" :{max:250,  hdr:"ALERT: SAFETY / LIMIT CONDITION",
    msg:"Risk: Flammable Gas\nPressure Exceeds Safe Limit!"},
  "C2H2":{max:17,   hdr:"ALERT: SAFETY / LIMIT CONDITION",
    msg:"Risk: Decomposition Hazard\nMax 17 bar for Acetylene!"},
  "H2S" :{max:20,   hdr:"ALERT: SAFETY / LIMIT CONDITION",
    msg:"Risk: Toxic Gas\nPressure Exceeds Safe Limit!"}
};
var GMAX = 2001;
var WG   = ["O2","H2","CH4","C2H2","H2S"];
var WMSG = {
  "O2" :"O\u2082 \u2014 Oxidizing gas. Use S14 O\u2082-cleaned equipment only.",
  "H2" :"H\u2082 \u2014 Highly flammable. S6 (\u2264200 bar) or S9 (\u2264400 bar) only.",
  "CH4":"CH\u2084 \u2014 Flammable. Max 250 bar. Use S16 double-contained spec.",
  "C2H2":"C\u2082H\u2082 \u2014 Unstable >17 bar. Decomposition risk!",
  "H2S":"H\u2082S \u2014 Toxic & flammable. Max 20 bar. PPE required."
};

function evalMode(gas, Pi) {
  if (Pi > GMAX)
    return {m:"D", hdr:"ALERT: SAFETY / LIMIT CONDITION",
      blocked:true, msg:"Risk: Pressure >2001 bar\nExceeds Maximum Limit!"};
  var lim = BLK[gas];
  if (lim && Pi > lim.max)
    return {m:"D", hdr:lim.hdr, blocked:true, msg:lim.msg};
  if (WG.indexOf(gas) !== -1)
    return {m:"W", hdr:"WARNING: PRESSURE LIMIT REACHED", blocked:false};
  return {m:"N", blocked:false};
}

function setMode(m, hdr) {
  var p  = document.getElementById("panel");
  p.className = m==="W"?"warn": m==="D"?"dngr":"";
  var ah = document.getElementById("alert-hdr");
  ah.style.display = m==="N" ? "none" : "flex";
  var at = document.getElementById("alert-txt");
  if (at && hdr) at.textContent = hdr;
}

// ── Correct toKey (parens removed, not replaced) ──────────────────────────────
function toKey(f){
  return f.replace(/ /g,"_").replace(/[()]/g,"").replace(/\//g,"").replace(/\u00b0/g,"deg");
}
function fv(l){var e=document.getElementById(toKey(l));return e?(parseFloat(e.value)||0):(DEF[l]||0);}
function getInletP(){
  var e=document.getElementById(toKey("Inlet Pressure (bar)"));
  return e?(parseFloat(e.value)||0):0;
}

// ── Tube lookup ────────────────────────────────────────────────────────────────
function lookupTube(gas,inP,diam){
  var all=TT.filter(function(r){
    return r.gas_codes.indexOf(gas)!==-1&&r.max_pressure>=inP&&r.id_mm>=diam;
  });
  if(!all.length)return null;
  var cap=all.filter(function(r){return r.max_pressure<=inP*5;});
  var c=cap.length?cap:all;
  c.sort(function(a,b){return(a.id_mm-b.id_mm)||(a.max_pressure-b.max_pressure);});
  return c[0];
}
function fmtTube(s){return s.replace(/X/g," \u00d7 ").replace(/\s{2,}/g," ").trim();}

function buildSpecCard(gas,inP,diam){
  var m=lookupTube(gas,inP,diam);
  if(!m)return"<div class='rbox' style='background:rgba(40,20,0,0.9);border-left:3px solid #ff8c00;color:#ffcc77;padding:8px 12px'><strong>\u26a0 No approved tube spec found</strong></div>";
  var badge="";
  if(gas==="O2")badge="<span class='badge b-o2'>O\u2082 \u2192 S14 only</span>";
  if(gas==="H2")badge="<span class='badge b-h2'>H\u2082 \u2192 S6/S9</span>";
  return"<div class='spec-card'>"
    +"<div class='sc-head'>\u25ba Recommended Tube Specification"+badge+"</div>"
    +"<div class='sc-row'>\u2713&nbsp;<strong>Tube Spec:&nbsp;&nbsp;&nbsp;&nbsp;</strong>"+m.spec+"<br>"
    +"\u2713&nbsp;<strong>Min. Tube Size:&nbsp;</strong>"+fmtTube(m.tube_od_w)+"</div>"
    +"<div class='sc-note'>ID "+m.id_mm.toFixed(3)+" mm&nbsp;\u2502&nbsp;Max rated: "+m.max_pressure+" bar</div>"
    +"</div>";
}

function calcConf(gas,inP,diam,spec){
  if(!spec)return 5;
  var dr=Math.min((spec.id_mm-diam)/spec.id_mm,1);
  var pr=Math.max((spec.max_pressure-inP)/spec.max_pressure,0);
  var sf=WG.indexOf(gas)!==-1?0.72:1.0;
  return Math.max(5,Math.min(100,Math.round((dr*0.35+pr*0.65)*sf*100)));
}

function showConf(pct){
  var w=document.getElementById("conf-wrap");
  w.className="show";
  setTimeout(function(){
    document.getElementById("conf-fill").style.width=pct+"%";
    document.getElementById("conf-pct").textContent=pct+"%";
  },80);
}

// ── Gas physics (Math.pow — no ** operator) ───────────────────────────────────
var GM={N2:.028013,O2:.031999,Ar:.039948,CO2:.04401,
        He:.0040026,H2:.002016,CH4:.01604,C2H2:.02604,
        FG1:.03881,FG2:.02671,Air:.02897};
var R=8.314,FR=0.02;
function rho(P,T,g){return(P*1e5*GM[g])/(R*(T+273.15));}
function rhoA(Pi,Po,T,g){return(rho(Pi,T,g)+rho(Po,T,g))/2;}
function calcDiameter(Pi,Po,T,L,Q,g){
  var dP=(Pi-Po)*1e5;if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  var Qs=Q/60000;
  return Math.pow((FR*L*8*rhoA(Pi,Po,T,g)*Qs*Qs)/(Math.pow(Math.PI,2)*dP),0.2)*1000;}
function calcFlow(Pi,Po,T,L,D,g){
  var dP=(Pi-Po)*1e5;if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  var Dm=D/1000;
  return Math.sqrt(dP*Math.pow(Math.PI,2)*Math.pow(Dm,5)/(8*FR*L*rhoA(Pi,Po,T,g)))*60000;}
function calcLength(Pi,Po,T,D,Q,g){
  var dP=(Pi-Po)*1e5;if(dP<=0)throw new Error("Inlet pressure must exceed outlet pressure.");
  var Dm=D/1000,Qs=Q/60000;
  return dP*Math.pow(Math.PI,2)*Math.pow(Dm,5)/(8*FR*rhoA(Pi,Po,T,g)*Qs*Qs);}
function calcOutlet(Pi,T,L,D,Q,g){
  var Dm=D/1000,Qs=Q/60000;
  function res(Po){return(Pi-Po)*1e5-(8*FR*L*rhoA(Pi,Po,T,g)*Qs*Qs)/(Math.pow(Math.PI,2)*Math.pow(Dm,5));}
  var lo=0,hi=Pi;
  for(var i=0;i<60;i++){if(Math.abs(hi-lo)<1e-4)break;var m=(lo+hi)/2;res(m)>0?lo=m:hi=m;}
  return Math.max((lo+hi)/2,0);}
function calcInlet(Po,T,L,D,Q,g){
  var lo=Po,hi=Po+10;
  while(hi<Po+2000){if(calcOutlet(hi,T,L,D,Q,g)>=Po)break;hi+=10;}
  for(var i=0;i<60;i++){var m=(lo+hi)/2,vm=calcOutlet(m,T,L,D,Q,g);
    if(Math.abs(vm-Po)<0.005)return m;vm<Po?lo=m:hi=m;}
  return(lo+hi)/2;}

// ── Fields ────────────────────────────────────────────────────────────────────
var FIELDS={
  diameter:["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  flow:    ["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  length:  ["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  inlet:   ["Temperature (\u00b0C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  outlet:  ["Temperature (\u00b0C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"]
};
var DEF={
  "Temperature (\u00b0C)":25,
  "Inlet Pressure (bar)":200,
  "Outlet Pressure (bar)":10,
  "Pipe Length (m)":100,
  "Pipe Diameter (mm)":10,
  "Flow Rate (LPM)":16
};

// ── Main calculation ───────────────────────────────────────────────────────────
function doCalc(){
  var gas=document.getElementById("gasSelect").value;
  var ct =document.getElementById("calcSelect").value;
  var th = evalMode(gas, ct==="inlet"?0:getInletP());
  setMode(th.m,th.hdr);

  var rN=document.getElementById("result-normal");
  var rC=document.getElementById("res-content");

  if(th.blocked&&ct!=="inlet"){
    rN.innerHTML="";
    rC.innerHTML="<div class='rbox blk'>\u26a0\ufe0f SAFETY BLOCKED<br>"
      +"<span style='font-weight:normal;font-size:12px'>"+th.msg.replace(/\n/g,"<br>")+"</span></div>";
    showConf(5);return;
  }
  try{
    var Tc=fv("Temperature (\u00b0C)"),Pi=fv("Inlet Pressure (bar)"),
        Po=fv("Outlet Pressure (bar)"),L=fv("Pipe Length (m)"),
        D=fv("Pipe Diameter (mm)"),Q=fv("Flow Rate (LPM)");
    var diam,inP,line;
    if(ct==="diameter"){var r=calcDiameter(Pi,Po,Tc,L,Q,gas);diam=r;inP=Pi;line="Required Diameter: <strong>"+r.toFixed(2)+" mm</strong>";}
    else if(ct==="flow")  {var r=calcFlow(Pi,Po,Tc,L,D,gas);  diam=D;inP=Pi;line="Maximum Flow Rate: <strong>"+r.toFixed(1)+" L/min</strong>";}
    else if(ct==="length"){var r=calcLength(Pi,Po,Tc,D,Q,gas);diam=D;inP=Pi;line="Maximum Pipe Length: <strong>"+r.toFixed(1)+" m</strong>";}
    else if(ct==="inlet") {
      var r=calcInlet(Po,Tc,L,D,Q,gas);diam=D;inP=r;
      var th2=evalMode(gas,r);setMode(th2.m,th2.hdr);
      if(th2.blocked){
        rN.innerHTML="";
        rC.innerHTML="<div style='font-size:13px;color:var(--tc2)'>Required Inlet: <strong>"+r.toFixed(2)+" bar</strong></div>"
          +"<div class='rbox blk' style='margin-top:6px'>\u26a0\ufe0f SAFETY BLOCKED<br>"
          +"<span style='font-weight:normal;font-size:12px'>"+th2.msg.replace(/\n/g,"<br>")+"</span></div>";
        showConf(10);return;
      }
      line="Required Inlet Pressure: <strong>"+r.toFixed(2)+" bar</strong>";
    }
    else{var r=calcOutlet(Pi,Tc,L,D,Q,gas);diam=D;inP=Pi;line="Estimated Outlet Pressure: <strong>"+r.toFixed(2)+" bar</strong>";}

    var spec=lookupTube(gas,inP,diam);
    var pct=calcConf(gas,inP,diam,spec);

    if(th.m==="N"){
      rN.innerHTML="<div class='rbox ok'><div class='res-title'>"+line+"</div></div>"
        +buildSpecCard(gas,inP,diam);
      rC.innerHTML="";
    } else {
      rN.innerHTML="";
      rC.innerHTML="<div style='font-size:13px;color:var(--tc2);line-height:1.85'>"+line+"<br>"
        +"Recommended Spec: <strong>"+(spec?spec.spec:"N/A")+"</strong><br>"
        +(spec?"Tube: <strong>"+fmtTube(spec.tube_od_w)+"</strong> (ID "+spec.id_mm.toFixed(2)+" mm)":"")
        +"</div>";
    }
    showConf(pct);
  }catch(e){
    document.getElementById("result-normal").innerHTML="<div class='rbox err'>\u26a0 "+e.message+"</div>";
    document.getElementById("res-content").innerHTML="";
  }
}

// ── Field builder ──────────────────────────────────────────────────────────────
function rebuildFields(){
  var gas=document.getElementById("gasSelect").value;
  var ct =document.getElementById("calcSelect").value;
  var fl =FIELDS[ct]||[];
  var html="";
  fl.forEach(function(f){
    var k=toKey(f);
    var extra=(f==="Inlet Pressure (bar)")?" oninput='onPressureInput()'":" ";
    html+="<div class='frow'><label>"+f+"</label>"
         +"<input type='number' id='"+k+"' value='"+(DEF[f]||0)+"' step='any'"+extra+"></div>";
  });
  document.getElementById("ifields").innerHTML=html;
  document.getElementById("result-normal").innerHTML="";
  document.getElementById("res-content").innerHTML="";
  document.getElementById("conf-wrap").className="";

  var Pi=ct==="inlet"?0:(DEF["Inlet Pressure (bar)"]||0);
  var th=evalMode(gas,Pi);
  setMode(th.m,th.hdr);
  if(th.m==="W"&&!th.blocked&&WMSG[gas]){
    document.getElementById("res-content").innerHTML=
      "<div style='font-size:12px;color:var(--tc2);opacity:.88'>"+WMSG[gas]+"</div>";
  }
}

function onPressureInput(){
  var gas=document.getElementById("gasSelect").value;
  var ct =document.getElementById("calcSelect").value;
  if(ct==="inlet")return;
  var th=evalMode(gas,getInletP());
  setMode(th.m,th.hdr);
}

rebuildFields();
</script>
</body>
</html>"""

page_html = PAGE.replace("%%HUD%%", HUD_TAG).replace("%%TUBE_TABLE%%", TUBE_TABLE_JSON)
components.html(page_html, height=950, scrolling=False)
