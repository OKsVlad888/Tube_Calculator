"""
High-Pressure Gas Flow Calculator  —  Upgraded UI Edition
==========================================================
Enhancements over base version:
  - Dynamic color theming: normal (teal) / warning (orange) / danger (red)
  - Safety alert header banners with animated glow
  - Recommendation Confidence progress bar
  - No-scroll compact layout (all content visible without scrolling)
  - Animated HUD background loaded from hud_background.html
  - Real-time theme update on gas/pressure change

Calculation logic and TUBE_TABLE are UNCHANGED from original.

Run:  streamlit run app.py
Req:  hud_background.html in the same directory
"""

import base64
import json
import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path

st.set_page_config(
    page_title="High-Pressure Gas Flow Calculator",
    page_icon="💨",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# Hide all Streamlit chrome; make iframe fill the viewport
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
#  TUBE-SPEC TABLE  —  Embedded directly (DO NOT MODIFY)
#
#  Schema per record:
#    gas_codes    : list of UI gas codes that this spec supports
#    spec         : tube spec name (S6, S9, S11, S12, S14, S16, B1, P6, C3)
#    id_mm        : inner diameter in mm
#    tube_od_w    : outer-diameter × wall-thickness label (e.g. '3/8"X0.049"')
#    max_pressure : maximum allowed working pressure [bar] for this row
# ══════════════════════════════════════════════════════════════════════════════

TUBE_TABLE = [
    # ── S6  SS-316  ≤200 bar ────────────────────────────────────────────────
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":4.9,    "tube_od_w":'1/4"X0.028"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":7.036,  "tube_od_w":'3/8"X0.049"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":9.398,  "tube_od_w":'1/2"X0.065"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":12.573, "tube_od_w":'5/8"X0.065"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":14.834, "tube_od_w":'3/4"X0.083"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":17.399, "tube_od_w":'7/8"X0.095"',"max_pressure":200.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S6","id_mm":19.863, "tube_od_w":'1"X0.109"',  "max_pressure":200.0},
    # ── S9  SS-316  ≤400 bar ────────────────────────────────────────────────
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":3.86,   "tube_od_w":'1/4"X0.049"',     "max_pressure":400.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":6.223,  "tube_od_w":'3/8"X0.065"',     "max_pressure":400.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":8.468,  "tube_od_w":'1/2"X0.083"',     "max_pressure":400.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":9.119,  "tube_od_w":'9/16"X0.10175"',  "max_pressure":400.0},
    {"gas_codes":["N2","Ar","CO2","He","CH4","C2H2","H2","FG1","FG2","Air"],"spec":"S9","id_mm":13.106, "tube_od_w":'3/4"X0.117"',     "max_pressure":400.0},
    # ── S11 SS-316  ≤2000 bar ───────────────────────────────────────────────
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S11","id_mm":2.108,  "tube_od_w":'1/4"X0.083"',"max_pressure":2000.0},
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S11","id_mm":3.1753, "tube_od_w":'3/8"X0.125"',"max_pressure":2000.0},
    # ── S12 SS-316  ≤1380 bar ───────────────────────────────────────────────
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S12","id_mm":2.6786, "tube_od_w":'1/4"X0.705"',   "max_pressure":1380.0},
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S12","id_mm":5.516,  "tube_od_w":'3/8"X0.0860"',  "max_pressure":1380.0},
    {"gas_codes":["N2","Ar","He","FG1","FG2","Air"],"spec":"S12","id_mm":7.925,  "tube_od_w":'9/16"X0.12525"',"max_pressure":1380.0},
    # ── S14 SS-316  (O2 approved, per-diameter pressure limits) ─────────────
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":4.572,  "tube_od_w":'1/4"X0.035"',"max_pressure":238.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":7.747,  "tube_od_w":'3/8"X0.035"',"max_pressure":170.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":10.21,  "tube_od_w":'1/2"X0.049"',"max_pressure":170.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":15.748, "tube_od_w":'3/4"X0.065"',"max_pressure":170.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":22.1,   "tube_od_w":'1"X0.065"',  "max_pressure":102.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":34.8,   "tube_od_w":'1.5"X0.065"',"max_pressure":68.0},
    {"gas_codes":["N2","O2","Ar","He","Air"],"spec":"S14","id_mm":47.5,   "tube_od_w":'2"X0.065"',  "max_pressure":61.0},
    # ── S16 SS-316 UHP Double-Contained  (per-diameter pressure limits) ─────
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":4.572,    "tube_od_w":'1/4"X0.035"', "max_pressure":204.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":7.747,    "tube_od_w":'3/8"X0.035"', "max_pressure":170.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":10.21,    "tube_od_w":'1/2"X0.049"', "max_pressure":170.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":13.38585, "tube_od_w":'5/8"X0.049"', "max_pressure":170.0},
    {"gas_codes":["CO2","CH4","C2H2","H2S"],"spec":"S16","id_mm":15.748,   "tube_od_w":'3/4"X0.065"', "max_pressure":170.0},
    # ── B1  Copper  ≤40 bar ─────────────────────────────────────────────────
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
    # ── P6  PPR  ≤10 bar ────────────────────────────────────────────────────
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
    # ── C3  Galvanized  ≤15 bar  (Air only) ─────────────────────────────────
    {"gas_codes":["Air"],"spec":"C3","id_mm":15.8,  "tube_od_w":'1/2" SCH40', "max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":15.59, "tube_od_w":'3/4" SCH40', "max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":26.24, "tube_od_w":'1" SCH40',   "max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":35.04, "tube_od_w":'1.25" SCH40',"max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":40.9,  "tube_od_w":'1.5" SCH40', "max_pressure":15.0},
    {"gas_codes":["Air"],"spec":"C3","id_mm":52.51, "tube_od_w":'2" SCH40',   "max_pressure":15.0},
]

# Serialize once for embedding in JS
TUBE_TABLE_JSON = json.dumps(TUBE_TABLE, ensure_ascii=False)


# ── Load optional HUD background HTML ────────────────────────────────────────
def load_hud() -> str:
    p = Path(__file__).resolve().parent / "hud_background.html"
    if not p.exists():
        return ""
    return "data:text/html;base64," + base64.b64encode(p.read_bytes()).decode()

hud_src = load_hud()
hud_tag = (
    f"<iframe id='hud' src='{hud_src}' sandbox='allow-scripts'></iframe>"
    if hud_src else
    "<div style='position:fixed;inset:0;background:#000d0d;z-index:1'></div>"
)


# ══════════════════════════════════════════════════════════════════════════════
#  FULL PAGE HTML — all CSS, HTML and JS inline
# ══════════════════════════════════════════════════════════════════════════════

PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
/* ── Reset & base ─────────────────────────────────────────────────────────── */
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
html,body{
  width:100%;height:100%;overflow:hidden;
  background:#000d0d;font-family:'Courier New',monospace;
  font-size:14px;color:#00e5cc;
}

/* ── CSS custom properties — theme tokens ────────────────────────────────── */
:root{
  --tc:#00e5cc;       /* primary theme color      */
  --tc2:#00ffee;      /* bright accent            */
  --tc-sub:rgba(0,229,204,0.42);    /* subtle border/bg   */
  --tc-glow:rgba(0,229,204,0.45);   /* glow color         */
  --tc-faint:rgba(0,229,204,0.08);  /* very faint fill    */
  --panel-bg:rgba(0,8,8,0.95);
  --panel-border:#00e5cc;
  --input-bg:rgba(0,14,12,0.98);
  --label-color:rgba(0,229,204,0.75);
  --btn-bg:rgba(0,229,204,0.07);
  --btn-hover:rgba(0,229,204,0.20);
  --glow-anim-start:rgba(0,229,204,0.28);
  --glow-anim-end:rgba(0,229,204,0.70);
  --alert-display:none;
  --conf-color:#00e5cc;
}

/* Orange theme — dangerous gas selected */
.theme-warning{
  --tc:#ffb800;
  --tc2:#ffe590;
  --tc-sub:rgba(255,184,0,0.42);
  --tc-glow:rgba(255,184,0,0.45);
  --tc-faint:rgba(255,184,0,0.06);
  --panel-bg:rgba(18,10,0,0.97);
  --panel-border:#ffb800;
  --input-bg:rgba(22,14,0,0.99);
  --label-color:rgba(255,200,80,0.80);
  --btn-bg:rgba(255,184,0,0.07);
  --btn-hover:rgba(255,184,0,0.22);
  --glow-anim-start:rgba(255,184,0,0.28);
  --glow-anim-end:rgba(255,184,0,0.72);
  --alert-display:flex;
  --conf-color:#ffb800;
}

/* Red theme — safety limit exceeded, calculation blocked */
.theme-danger{
  --tc:#ff3333;
  --tc2:#ff8888;
  --tc-sub:rgba(255,60,60,0.42);
  --tc-glow:rgba(255,60,60,0.45);
  --tc-faint:rgba(255,60,60,0.06);
  --panel-bg:rgba(18,0,0,0.98);
  --panel-border:#ff3333;
  --input-bg:rgba(22,2,2,0.99);
  --label-color:rgba(255,140,140,0.80);
  --btn-bg:rgba(255,60,60,0.07);
  --btn-hover:rgba(255,60,60,0.22);
  --glow-anim-start:rgba(255,60,60,0.28);
  --glow-anim-end:rgba(255,60,60,0.72);
  --alert-display:flex;
  --conf-color:#ff3333;
}

/* ── HUD iframe ───────────────────────────────────────────────────────────── */
#hud{
  position:fixed;inset:0;width:100%;height:100%;
  border:none;pointer-events:none;z-index:1;
}

/* ── Outer wrapper — centers the panel ───────────────────────────────────── */
#wrapper{
  position:fixed;top:0;left:0;width:100%;height:100%;
  display:flex;align-items:center;justify-content:center;
  z-index:500;pointer-events:none;
  padding:12px;
}

/* ── Main panel ───────────────────────────────────────────────────────────── */
#panel{
  pointer-events:all;
  width:min(98%,560px);
  max-height:96vh;
  overflow:hidden;       /* NO internal scrollbar */
  background:var(--panel-bg);
  border:1.5px solid var(--panel-border);
  border-radius:4px;
  padding:10px 20px 14px;
  position:relative;
  animation:panelGlow 2.5s ease-in-out infinite alternate;
  transition:background 0.4s, border-color 0.4s;
}

/* Panel glow animation — uses theme tokens */
@keyframes panelGlow{
  from{box-shadow:0 0 8px var(--glow-anim-start),inset 0 0 20px var(--tc-faint)}
  to  {box-shadow:0 0 32px var(--glow-anim-end), inset 0 0 50px var(--tc-faint)}
}

/* CRT scanline overlay on panel */
#panel::after{
  content:"";position:absolute;inset:0;pointer-events:none;border-radius:4px;z-index:0;
  background:repeating-linear-gradient(
    to bottom,transparent 0,transparent 3px,
    rgba(0,229,204,0.008) 3px,rgba(0,229,204,0.008) 4px);
}

/* ── Alert header (hidden in normal mode) ────────────────────────────────── */
#alert-header{
  display:var(--alert-display);
  align-items:center;justify-content:center;gap:10px;
  padding:6px 12px;margin-bottom:8px;
  background:var(--tc-faint);
  border:1px solid var(--tc-sub);
  border-radius:3px;
  position:relative;z-index:1;
  animation:alertPulse 1.2s ease-in-out infinite alternate;
  transition:display 0.3s;
}
@keyframes alertPulse{
  from{box-shadow:0 0 4px var(--glow-anim-start)}
  to  {box-shadow:0 0 18px var(--glow-anim-end)}
}
#alert-header .alert-icon{font-size:15px}
#alert-header .alert-title{
  color:var(--tc);font-size:12px;
  font-weight:bold;letter-spacing:2.5px;text-transform:uppercase;
  text-shadow:0 0 10px var(--tc-glow);
}

/* ── Panel title ─────────────────────────────────────────────────────────── */
#panel h1{
  position:relative;z-index:1;
  color:var(--tc2);text-shadow:0 0 14px var(--tc-glow);
  font-size:13px;letter-spacing:3.5px;text-transform:uppercase;text-align:center;
  margin-bottom:10px;line-height:1.45;
  transition:color 0.4s, text-shadow 0.4s;
}
#panel h1 .tri{color:var(--tc);margin-right:7px}

/* ── Body container ───────────────────────────────────────────────────────── */
#body{position:relative;z-index:1;display:flex;flex-direction:column;gap:7px}

/* ── Select rows ─────────────────────────────────────────────────────────── */
.sr{display:flex;flex-direction:column;gap:3px}
.sr label{
  font-size:12px;letter-spacing:1.5px;text-transform:uppercase;
  color:var(--label-color);transition:color 0.4s;
}
.sr select{
  background:var(--input-bg);color:var(--tc2);
  border:1px solid var(--tc-sub);border-radius:2px;
  padding:5px 28px 5px 8px;
  font-family:'Courier New',monospace;font-size:13px;
  width:100%;cursor:pointer;outline:none;
  -webkit-appearance:none;appearance:none;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%2300e5cc'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 10px center;
  transition:background-color 0.4s,border-color 0.4s,color 0.4s;
}
.sr select:focus{border-color:var(--tc);box-shadow:0 0 7px var(--tc-glow)}
.sr select option{background:#000d0d;color:#00e5cc}

/* ── Divider ─────────────────────────────────────────────────────────────── */
.hdiv{border:none;border-top:1px solid var(--tc-sub);margin:2px 0;transition:border-color 0.4s}

/* ── Input fields ────────────────────────────────────────────────────────── */
#input-fields{display:flex;flex-direction:column;gap:5px}
.frow{display:flex;flex-direction:column;gap:2px}
.frow label{font-size:12px;color:var(--label-color);letter-spacing:0.2px;transition:color 0.4s}
.frow input{
  width:100%;background:var(--input-bg);color:var(--tc2);
  border:1px solid var(--tc-sub);border-radius:2px;
  padding:5px 8px;font-family:'Courier New',monospace;font-size:13px;outline:none;
  transition:background-color 0.4s,border-color 0.4s,color 0.4s;
}
.frow input:focus{border-color:var(--tc);box-shadow:0 0 5px var(--tc-glow)}
.frow input::-webkit-inner-spin-button,
.frow input::-webkit-outer-spin-button{-webkit-appearance:none;margin:0}
.frow input[type=number]{-moz-appearance:textfield}

/* ── Caption ─────────────────────────────────────────────────────────────── */
.cap{font-size:11px;color:rgba(0,229,204,0.38);letter-spacing:0.8px}

/* ── Calculate button ────────────────────────────────────────────────────── */
#cbtn{
  width:100%;background:var(--btn-bg);color:var(--tc2);
  border:1.5px solid var(--tc);border-radius:2px;padding:7px 0;
  font-family:'Courier New',monospace;font-size:13px;
  letter-spacing:4px;cursor:pointer;text-transform:uppercase;
  transition:background 0.2s,box-shadow 0.2s,border-color 0.4s,color 0.4s;
}
#cbtn:hover{background:var(--btn-hover);box-shadow:0 0 16px var(--tc-glow)}
#cbtn:active{background:rgba(255,255,255,0.08)}

/* ── Result boxes ────────────────────────────────────────────────────────── */
.rbox{padding:8px 12px;border-radius:2px;font-size:13px;line-height:1.7;word-break:break-word}
.rbox.ok{
  background:rgba(0,56,28,0.94);border:1px solid rgba(0,200,100,0.40);
  border-left:3px solid #00e87a;color:#00ffcc;
}
.rbox.ok .res-title{font-size:14px;font-weight:bold}
.rbox.warn{
  background:rgba(50,30,0,0.96);border:1px solid rgba(255,180,40,0.35);
  border-left:3px solid #ffb800;color:#ffe590;
}
.rbox.danger{
  background:rgba(50,0,0,0.97);border:1px solid rgba(255,60,60,0.35);
  border-left:3px solid #ff3333;color:#ff8080;
}
.rbox.err{
  background:rgba(38,0,0,0.96);border:1px solid rgba(255,60,60,0.32);
  border-left:3px solid #ff4444;color:#ff9090;
}

/* Safety blocked result box */
.rbox.blocked{
  background:rgba(50,0,0,0.97);border:1px solid rgba(255,60,60,0.40);
  border-left:3px solid #ff3333;color:#ff8080;font-weight:bold;
  animation:blockPulse 1.5s ease-in-out infinite alternate;
}
@keyframes blockPulse{
  from{box-shadow:0 0 4px rgba(255,60,60,0.2)}
  to  {box-shadow:0 0 14px rgba(255,60,60,0.5)}
}

.rbox .he{direction:rtl;font-size:12px;opacity:.80;margin-top:3px}

/* ── Spec card ───────────────────────────────────────────────────────────── */
.spec-card{
  margin-top:6px;padding:8px 12px;border-radius:2px;
  background:rgba(0,28,22,0.98);
  border:1px solid var(--tc-sub);
  border-left:3px solid var(--tc);
  transition:border-color 0.4s;
}
.spec-card .sc-head{
  color:rgba(0,229,204,0.52);font-size:10px;
  letter-spacing:2px;text-transform:uppercase;margin-bottom:5px;
}
.spec-card .sc-row{color:#00ffee;font-size:13px;line-height:1.9}
.spec-card .sc-note{color:rgba(0,229,204,0.45);font-size:11px;margin-top:2px}
.badge{
  display:inline-block;border-radius:2px;font-size:9px;
  letter-spacing:1.5px;padding:1px 6px;margin-left:6px;
  vertical-align:middle;text-transform:uppercase;
}
.b-o2{background:rgba(255,100,0,0.15);border:1px solid rgba(255,120,0,0.4);color:#ffaa44}
.b-h2{background:rgba(0,180,255,0.12);border:1px solid rgba(0,200,255,0.35);color:#66ddff}

/* ── Confidence bar ──────────────────────────────────────────────────────── */
#conf-container{
  margin-top:6px;padding:6px 10px;
  background:rgba(0,0,0,0.3);
  border:1px solid var(--tc-sub);border-radius:2px;
  display:none;          /* shown after first calculation */
  flex-direction:column;gap:4px;
  transition:border-color 0.4s;
}
#conf-container.show{display:flex}
.conf-label{
  font-size:9px;letter-spacing:2px;text-transform:uppercase;
  color:var(--label-color);
}
.conf-track{
  position:relative;height:10px;
  background:rgba(0,0,0,0.4);border:1px solid var(--tc-sub);border-radius:2px;
  overflow:hidden;
}
.conf-fill{
  height:100%;background:var(--tc);
  transition:width 0.6s ease,background 0.4s;
  position:relative;
}
.conf-fill::after{
  content:'';position:absolute;right:0;top:0;width:3px;height:100%;
  background:var(--tc2);opacity:0.8;
  box-shadow:0 0 6px var(--tc-glow);
}
.conf-pct{
  font-size:11px;color:var(--tc);font-weight:bold;letter-spacing:1px;
  text-align:right;transition:color 0.4s;
}
</style>
</head>
<body>

%%HUD%%

<div id="wrapper">
  <div id="panel">

    <!-- Safety alert header — shown in warning/danger modes -->
    <div id="alert-header">
      <span class="alert-icon">&#9888;</span>
      <span class="alert-title" id="alert-title">WARNING</span>
    </div>

    <!-- Calculator title -->
    <h1><span class="tri">&#9651;</span>HIGH-PRESSURE GAS FLOW CALCULATOR</h1>

    <div id="body">

      <!-- Gas type selector -->
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

      <!-- Calculation type selector -->
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

      <!-- Dynamic input fields -->
      <div id="input-fields"></div>

      <p class="cap">Friction factor (f) = 0.02 &nbsp;&middot;&nbsp; fixed constant</p>

      <button id="cbtn" onclick="doCalc()">&#9658;&nbsp;&nbsp;Calculate</button>

      <!-- Results area -->
      <div id="result-area"></div>

      <!-- Recommendation confidence bar -->
      <div id="conf-container">
        <div class="conf-label">Recommendation Confidence</div>
        <div class="conf-track">
          <div class="conf-fill" id="conf-fill" style="width:0%"></div>
        </div>
        <div class="conf-pct" id="conf-pct">0%</div>
      </div>

    </div><!-- /#body -->
  </div><!-- /#panel -->
</div><!-- /#wrapper -->

<script>
// ═══════════════════════════════════════════════════════════════════════════
//  TUBE_TABLE  —  injected from Python (DO NOT MODIFY)
// ═══════════════════════════════════════════════════════════════════════════
var TUBE_TABLE = %%TUBE_TABLE%%;

// ═══════════════════════════════════════════════════════════════════════════
//  SAFETY / THEME CONFIGURATION
// ═══════════════════════════════════════════════════════════════════════════

// Gases that always trigger orange warning panel
var DANGEROUS_GASES = ["O2","H2","CH4","C2H2","H2S"];

// Gases with a specific inlet-pressure limit that triggers red / blocked
var SAFETY_LIMITS = {
  "CO2":  { max_bar:69.9,
    title:"ALERT: SAFETY / LIMIT CONDITION",
    msg:"Risk: Liquid Phase CO\u2082\nPressure Exceeds Limit!",
    he:"\u05dc\u05d7\u05e5 \u05db\u05e0\u05d9\u05e1\u05d4 \u05d7\u05d5\u05e8\u05d2 \u2014 CO\u2082 \u05e2\u05d5\u05d1\u05e8 \u05dc\u05e4\u05d0\u05d6\u05d4 \u05e0\u05d5\u05d6\u05dc\u05d9\u05ea"
  },
  "CH4":  { max_bar:250,
    title:"ALERT: SAFETY / LIMIT CONDITION",
    msg:"Risk: Flammable Gas\nPressure Exceeds Safe Limit!",
    he:"\u05dc\u05d7\u05e5 \u05db\u05e0\u05d9\u05e1\u05d4 \u05d7\u05d5\u05e8\u05d2 \u2014 \u05de\u05ea\u05d0\u05df \u05d1\u05dc\u05d7\u05e5 \u05de\u05e7\u05e1\u05d9\u05de\u05dc\u05d9"
  },
  "C2H2": { max_bar:17,
    title:"ALERT: SAFETY / LIMIT CONDITION",
    msg:"Risk: Decomposition Hazard\nMax 17 bar for Acetylene!",
    he:"\u05d1\u05d8\u05d9\u05d7\u05d5\u05ea \u2014 \u05dc\u05d7\u05e5 \u05de\u05e7\u05e1\u05d9\u05de\u05dc\u05d9 \u05dc\u05d0\u05e6\u05d8\u05d9\u05dc\u05df: 17 \u05d1\u05e8"
  },
  "H2S":  { max_bar:20,
    title:"ALERT: SAFETY / LIMIT CONDITION",
    msg:"Risk: Toxic Gas\nPressure Exceeds Safe Limit!",
    he:"\u05d1\u05d8\u05d9\u05d7\u05d5\u05ea \u2014 \u05dc\u05d7\u05e5 \u05de\u05e7\u05e1\u05d9\u05de\u05dc\u05d9 \u05dc-H\u2082S: 20 \u05d1\u05e8"
  }
};

// Warning messages for orange (dangerous gas) state
var DANGER_WARNINGS = {
  "O2":  "O\u2082 \u2014 Oxidizing gas. Use O\u2082-cleaned fittings and S14 spec only.",
  "H2":  "H\u2082 \u2014 Highly flammable. Ensure leak-free connections and S6/S9 spec.",
  "CH4": "CH\u2084 \u2014 Flammable gas. Ensure adequate ventilation.",
  "C2H2":"C\u2082H\u2082 \u2014 Unstable >17 bar. Keep pressure below decomposition limit.",
  "H2S": "H\u2082S \u2014 Highly toxic. Mandatory gas detection and PPE required."
};

// ═══════════════════════════════════════════════════════════════════════════
//  THEME ENGINE
// ═══════════════════════════════════════════════════════════════════════════

/**
 * Apply visual theme to the panel.
 * @param {string} mode  'normal' | 'warning' | 'danger'
 * @param {string} title  Alert header title text (used in warning/danger)
 */
function setTheme(mode, title) {
  var panel = document.getElementById("panel");
  panel.className = mode === "normal" ? "" : ("theme-" + mode);

  var header = document.getElementById("alert-header");
  var titleEl = document.getElementById("alert-title");

  if (mode === "normal") {
    header.style.display = "none";
  } else {
    header.style.display = "flex";
    if (titleEl) titleEl.textContent = title || "WARNING";
  }
}

/** Read current inlet pressure from the input field (if visible). */
function getCurrentInletP() {
  // toKey("Inlet Pressure (bar)") → "Inlet_Pressure_bar_"
  var el = document.getElementById("Inlet_Pressure_bar_");
  return el ? (parseFloat(el.value) || 0) : 0;
}

/**
 * Determine which theme should apply given gas + inlet pressure.
 * Returns {mode, title, blocked, blockInfo}
 */
function evalTheme(gas, inletP) {
  // Safety block (red) takes priority
  var lim = SAFETY_LIMITS[gas];
  if (lim && inletP > lim.max_bar) {
    return { mode:"danger", title:lim.title, blocked:true, blockInfo:lim };
  }
  // Dangerous gas warning (orange)
  if (DANGEROUS_GASES.indexOf(gas) !== -1) {
    return {
      mode:"warning",
      title:"WARNING: PRESSURE LIMIT REACHED",
      blocked:false
    };
  }
  return { mode:"normal", blocked:false };
}

/** Called on gas select change — updates theme but does not run calc. */
function onGasChange() {
  rebuildFields();
  var gas = document.getElementById("gasSelect").value;
  var P   = getCurrentInletP();
  var t   = evalTheme(gas, P);
  setTheme(t.mode, t.title);
  // Show inline gas warning in result area (informational only)
  if (t.mode === "warning" && DANGER_WARNINGS[gas]) {
    document.getElementById("result-area").innerHTML =
      "<div class='rbox warn'><strong>&#9888;&nbsp;" + gas + " — Hazardous Gas</strong><br>" +
      DANGER_WARNINGS[gas] + "</div>";
  } else if (t.mode === "normal") {
    document.getElementById("result-area").innerHTML = "";
  }
  hideConfBar();
}

/** Called on inlet pressure change — live theme update. */
function onPressureChange() {
  var gas = document.getElementById("gasSelect").value;
  var P   = getCurrentInletP();
  var t   = evalTheme(gas, P);
  setTheme(t.mode, t.title);
}

// ═══════════════════════════════════════════════════════════════════════════
//  CONFIDENCE BAR
// ═══════════════════════════════════════════════════════════════════════════

/** Compute 0–100 recommendation confidence score. */
function calcConfidence(gas, inletP, requiredD, tubeSpec) {
  if (!tubeSpec) return 5;

  // Diameter margin: how much bigger the tube ID is vs required (0–1 capped)
  var dRatio  = Math.min((tubeSpec.id_mm - requiredD) / tubeSpec.id_mm, 1);
  // Pressure margin: headroom below rated max (0–1)
  var pRatio  = Math.max((tubeSpec.max_pressure - inletP) / tubeSpec.max_pressure, 0);
  // Safety penalty for dangerous gases
  var safetyF = (DANGEROUS_GASES.indexOf(gas) !== -1) ? 0.72 : 1.0;

  var score   = (dRatio * 0.35 + pRatio * 0.65) * safetyF;
  return Math.max(5, Math.min(100, Math.round(score * 100)));
}

function showConfBar(pct) {
  var container = document.getElementById("conf-container");
  var fill      = document.getElementById("conf-fill");
  var label     = document.getElementById("conf-pct");
  container.className = "show";
  // Use setTimeout to allow CSS transition
  setTimeout(function() {
    fill.style.width = pct + "%";
    label.textContent = pct + "%";
  }, 60);
}

function hideConfBar() {
  document.getElementById("conf-container").className = "";
}

// ═══════════════════════════════════════════════════════════════════════════
//  TUBE TABLE LOOKUP  —  unchanged from original
// ═══════════════════════════════════════════════════════════════════════════
var OVER_ENG_CAP = 5;

function lookupTubeSpec(gas, inletP, diam) {
  var all = TUBE_TABLE.filter(function(r) {
    return r.gas_codes.indexOf(gas) !== -1 &&
           r.max_pressure >= inletP &&
           r.id_mm >= diam;
  });
  if (!all.length) return null;

  var capped = all.filter(function(r) { return r.max_pressure <= inletP * OVER_ENG_CAP; });
  var c = capped.length ? capped : all;

  c.sort(function(a, b) {
    return (a.id_mm - b.id_mm) || (a.max_pressure - b.max_pressure);
  });
  return c[0];
}

function fmtTube(s) {
  return s.replace(/X/g," \u00d7 ").replace(/\s{2,}/g," ").trim();
}

function buildSpecCard(gas, inletP, diam) {
  var m = lookupTubeSpec(gas, inletP, diam);
  if (!m) {
    return "<div class='rbox warn' style='margin-top:6px'>"
      + "<strong>&#9888;&nbsp;No approved tube specification found</strong><br>"
      + "<span class='he'>\u05d0\u05d9\u05df \u05de\u05e4\u05e8\u05d8 \u05e6\u05e0\u05e8\u05ea \u05de\u05d0\u05d5\u05e9\u05e8</span></div>";
  }
  var badge = "";
  if (gas === "O2") badge = "<span class='badge b-o2'>O\u2082 \u2192 S14 only</span>";
  if (gas === "H2") badge = "<span class='badge b-h2'>H\u2082 \u2192 S6 / S9</span>";
  return "<div class='spec-card'>"
    + "<div class='sc-head'>&#9654; Recommended Tube Specification" + badge + "</div>"
    + "<div class='sc-row'>"
    + "&#10003;&nbsp;<strong>Tube Spec:&nbsp;&nbsp;&nbsp;&nbsp;</strong>" + m.spec + "<br>"
    + "&#10003;&nbsp;<strong>Min. Tube Size:&nbsp;</strong>" + fmtTube(m.tube_od_w)
    + "</div>"
    + "<div class='sc-note'>"
    + "ID " + m.id_mm.toFixed(3) + " mm"
    + "&nbsp;&nbsp;\u2502&nbsp;&nbsp;"
    + "Max pressure rated: " + m.max_pressure + " bar"
    + "</div></div>";
}

function renderResult(line, diam, inletP, gas) {
  var spec = lookupTubeSpec(gas, inletP, diam);
  var pct  = calcConfidence(gas, inletP, diam, spec);
  showConfBar(pct);
  return "<div class='rbox ok'><div class='res-title'>" + line + "</div></div>"
       + buildSpecCard(gas, inletP, diam);
}

// ═══════════════════════════════════════════════════════════════════════════
//  GAS PHYSICS  —  UNCHANGED from original
// ═══════════════════════════════════════════════════════════════════════════
var GAS_M = {
  N2:.028013, O2:.031999, Ar:.039948, CO2:.04401,
  He:.0040026, H2:.002016, CH4:.01604, C2H2:.02604,
  FG1:.03881, FG2:.02671, Air:.02897
};
var R = 8.314, FR = 0.02;

function rho(P, T, g)    { return (P*1e5*GAS_M[g]) / (R*(T+273.15)); }
function rhoA(Pi,Po,T,g) { return (rho(Pi,T,g)+rho(Po,T,g))/2; }

function calcDiameter(Pi,Po,T,L,Q,g){
  var dP=(Pi-Po)*1e5; if(dP<=0) throw new Error("Inlet pressure must exceed outlet pressure.");
  return Math.pow((FR*L*8*rhoA(Pi,Po,T,g)*(Q/60000)**2)/(Math.PI**2*dP),0.2)*1000;
}
function calcFlow(Pi,Po,T,L,D,g){
  var dP=(Pi-Po)*1e5; if(dP<=0) throw new Error("Inlet pressure must exceed outlet pressure.");
  return Math.sqrt(dP*Math.PI**2*(D/1000)**5/(8*FR*L*rhoA(Pi,Po,T,g)))*60000;
}
function calcLength(Pi,Po,T,D,Q,g){
  var dP=(Pi-Po)*1e5; if(dP<=0) throw new Error("Inlet pressure must exceed outlet pressure.");
  return dP*Math.PI**2*(D/1000)**5/(8*FR*rhoA(Pi,Po,T,g)*(Q/60000)**2);
}
function calcOutlet(Pi,T,L,D,Q,g){
  var Dm=D/1000, Qs=Q/60000;
  function res(Po){return(Pi-Po)*1e5-(8*FR*L*rhoA(Pi,Po,T,g)*Qs**2)/(Math.PI**2*Dm**5)}
  var lo=0, hi=Pi;
  for(var i=0;i<60;i++){if(Math.abs(hi-lo)<1e-4)break;var m=(lo+hi)/2;res(m)>0?lo=m:hi=m}
  return Math.max((lo+hi)/2,0);
}
function calcInlet(Po,T,L,D,Q,g){
  var lo=Po, hi=Po+10;
  while(hi<Po+2000){if(calcOutlet(hi,T,L,D,Q,g)>=Po)break;hi+=10}
  for(var i=0;i<60;i++){
    var m=(lo+hi)/2, vm=calcOutlet(m,T,L,D,Q,g);
    if(Math.abs(vm-Po)<0.005)return m;
    vm<Po?lo=m:hi=m;
  }
  return(lo+hi)/2;
}

// ═══════════════════════════════════════════════════════════════════════════
//  FIELD CONFIGURATION  —  unchanged
// ═══════════════════════════════════════════════════════════════════════════
var FIELDS = {
  diameter: ["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  flow:     ["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  length:   ["Temperature (\u00b0C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  inlet:    ["Temperature (\u00b0C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  outlet:   ["Temperature (\u00b0C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"]
};
var DEF = {
  "Temperature (\u00b0C)":25,
  "Inlet Pressure (bar)":100,
  "Outlet Pressure (bar)":10,
  "Pipe Length (m)":10,
  "Pipe Diameter (mm)":10,
  "Flow Rate (LPM)":100
};

// Convert field label to a safe HTML id
function toKey(f) {
  return f.replace(/ /g,"_").replace(/[()]/g,"_").replace(/\//g,"").replace(/\u00b0/g,"deg");
}

// Read a field value by label
function fv(l) {
  var e = document.getElementById(toKey(l));
  return e ? (parseFloat(e.value)||0) : (DEF[l]||0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  MAIN CALCULATION  —  updated with safety-block and theming
// ═══════════════════════════════════════════════════════════════════════════
function doCalc() {
  var gas  = document.getElementById("gasSelect").value;
  var ct   = document.getElementById("calcSelect").value;
  var area = document.getElementById("result-area");

  try {
    var Tc = fv("Temperature (\u00b0C)");
    var Pi = fv("Inlet Pressure (bar)");
    var Po = fv("Outlet Pressure (bar)");
    var L  = fv("Pipe Length (m)");
    var D  = fv("Pipe Diameter (mm)");
    var Q  = fv("Flow Rate (LPM)");

    // Determine inlet for theme check (unknown for 'inlet' calc type)
    var checkP = (ct === "inlet") ? 0 : Pi;
    var theme  = evalTheme(gas, checkP);
    setTheme(theme.mode, theme.title);

    // ── Safety block: refuse calculation and show blocked result ──────────
    if (theme.blocked && ct !== "inlet") {
      var bi = theme.blockInfo;
      area.innerHTML =
        "<div class='rbox blocked'>"
        + "<strong>&#9888; SAFETY BLOCKED</strong><br>"
        + bi.msg.replace(/\n/g,"<br>")
        + "<div class='he'>" + bi.he + "</div>"
        + "</div>";
      showConfBar(5);
      return;
    }

    // ── Normal calculation ────────────────────────────────────────────────
    var html = "";

    if (ct === "diameter") {
      var r = calcDiameter(Pi,Po,Tc,L,Q,gas);
      html = renderResult("Required Diameter: <strong>" + r.toFixed(2) + " mm</strong>", r, Pi, gas);

    } else if (ct === "flow") {
      var r = calcFlow(Pi,Po,Tc,L,D,gas);
      html = renderResult("Maximum Flow Rate: <strong>" + r.toFixed(1) + " L/min</strong>", D, Pi, gas);

    } else if (ct === "length") {
      var r = calcLength(Pi,Po,Tc,D,Q,gas);
      html = renderResult("Maximum Pipe Length: <strong>" + r.toFixed(1) + " m</strong>", D, Pi, gas);

    } else if (ct === "inlet") {
      var r = calcInlet(Po,Tc,L,D,Q,gas);
      // Re-check safety with calculated inlet
      var t2 = evalTheme(gas, r);
      setTheme(t2.mode, t2.title);
      if (t2.blocked) {
        var bi = t2.blockInfo;
        html = "<div class='rbox ok'><div class='res-title'>Required Inlet Pressure: <strong>"
             + r.toFixed(2) + " bar</strong></div></div>"
             + "<div class='rbox blocked' style='margin-top:6px'>"
             + "<strong>&#9888; Calculated pressure exceeds safe limit for " + gas + "</strong><br>"
             + bi.msg.replace(/\n/g,"<br>")
             + "<div class='he'>" + bi.he + "</div>"
             + "</div>";
        showConfBar(10);
      } else {
        html = renderResult("Required Inlet Pressure: <strong>" + r.toFixed(2) + " bar</strong>", D, r, gas);
      }

    } else if (ct === "outlet") {
      var r = calcOutlet(Pi,Tc,L,D,Q,gas);
      html = renderResult("Estimated Outlet Pressure: <strong>" + r.toFixed(2) + " bar</strong>", D, Pi, gas);
    }

    area.innerHTML = html;

  } catch(e) {
    area.innerHTML = "<div class='rbox err'>&#9888; " + e.message + "</div>";
    hideConfBar();
  }
}

// ═══════════════════════════════════════════════════════════════════════════
//  FIELD BUILDER  —  adds oninput handler to inlet pressure for live theme
// ═══════════════════════════════════════════════════════════════════════════
function rebuildFields() {
  var ct = document.getElementById("calcSelect").value;
  var fl = FIELDS[ct] || [];
  var html = "";
  fl.forEach(function(f) {
    var k = toKey(f);
    var extra = (f === "Inlet Pressure (bar)") ? " oninput='onPressureChange()'" : "";
    html += "<div class='frow'><label>" + f + "</label>"
          + "<input type='number' id='" + k + "' value='" + (DEF[f]||0) + "' step='any'" + extra + "></div>";
  });
  document.getElementById("input-fields").innerHTML = html;
  document.getElementById("result-area").innerHTML = "";
  hideConfBar();
}

// ── Initialise on load ───────────────────────────────────────────────────
rebuildFields();
</script>
</body>
</html>"""

page_html = PAGE.replace("%%HUD%%", hud_tag).replace("%%TUBE_TABLE%%", TUBE_TABLE_JSON)
components.html(page_html, height=950, scrolling=False)
