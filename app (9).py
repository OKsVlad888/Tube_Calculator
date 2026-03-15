"""
High-Pressure Gas Flow Calculator
• hud_background.html  → full-screen animated background
• Title: centered between the two panels (top strip)
• Left panel : Gas type | Calc type | first 4 input fields  (fits content, scrollable)
• Right panel: last input field | CALCULATE | result box | spec table  (fits content, scrollable)
• HUD graphics visible: top strip, bottom strip, far-left, far-right, centre gap
• Run: streamlit run app.py   (hud_background.html in same folder)
"""

import math, base64
import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path

st.set_page_config(page_title="Gas Flow Calculator", page_icon="💨",
                   layout="wide", initial_sidebar_state="collapsed")

# ── Gas / calc data ────────────────────────────────────────────────────────
GAS_DATA = {
    "N2":0.028013,"O2":0.031999,"Ar":0.039948,"CO2":0.04401,
    "He":0.0040026,"H2":0.002016,"CH4":0.01604,"C2H2":0.02604,
    "Forming Gas 1":0.03881,"Forming Gas 2":0.02671,"Air":0.02897,
}
GAS_LABEL = {
    "N2":"N2 (Nitrogen)","O2":"O2 (Oxygen)","Ar":"Ar (Argon)",
    "CO2":"CO2 (Carbon Dioxide)","He":"He (Helium)","H2":"H2 (Hydrogen)",
    "CH4":"CH4 (Methane)","C2H2":"C2H2 (Acetylene)",
    "Forming Gas 1":"H2 3% + Ar 97%  (FG1)",
    "Forming Gas 2":"H2 5% + N2 95%  (FG2)",
    "Air":"Air (Dry Air)",
}
CALC_FIELDS = {
    "Pipe Diameter (mm)":    ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
    "Flow Rate (LPM)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
    "Pipe Length (m)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
    "Inlet Pressure (bar)":  ["Temperature (°C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
    "Outlet Pressure (bar)": ["Temperature (°C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
}
DEFAULTS = {
    "Temperature (°C)":25.0,"Inlet Pressure (bar)":100.0,
    "Outlet Pressure (bar)":10.0,"Pipe Length (m)":10.0,
    "Pipe Diameter (mm)":10.0,"Flow Rate (LPM)":100.0,
}

# ── Physics ────────────────────────────────────────────────────────────────
def _rho(P,T,g):    return (P*GAS_DATA[g])/(8.314*T)
def _ra(Pi,Po,T,g): return (_rho(Pi*1e5,T,g)+_rho(Po*1e5,T,g))/2
FR=0.02

def calc_diameter(Pi,Po,Tc,L,Q,g):
    T=Tc+273.15; Qs=Q/6e4; dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return ((FR*L*8*_ra(Pi,Po,T,g)*Qs**2)/(math.pi**2*dP))**0.2*1000

def calc_flow(Pi,Po,Tc,L,D,g):
    T=Tc+273.15; Dm=D/1000; dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return math.sqrt((dP*math.pi**2*Dm**5)/(8*FR*L*_ra(Pi,Po,T,g)))*6e4

def calc_length(Pi,Po,Tc,D,Q,g):
    T=Tc+273.15; Dm=D/1000; Qs=Q/6e4; dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return (dP*math.pi**2*Dm**5)/(8*FR*_ra(Pi,Po,T,g)*Qs**2)

def calc_outlet(Pi,Tc,L,D,Q,g):
    """Bisection solver – guaranteed convergence, no oscillation."""
    T=Tc+273.15; Dm=D/1000; Qs=Q/6e4
    def residual(Po):
        ra=_ra(Pi,Po,T,g)
        dP_calc=(8*FR*L*ra*Qs**2)/(math.pi**2*Dm**5)
        return (Pi-Po)*1e5 - dP_calc   # zero when Po is correct
    lo,hi=0.0,Pi
    for _ in range(60):
        if abs(hi-lo)<1e-4: break
        mid=(lo+hi)/2
        lo,hi=(mid,hi) if residual(mid)>0 else (lo,mid)
    return max((lo+hi)/2, 0.0)

def calc_inlet(Po,Tc,L,D,Q,g):
    # outlet pressure increases with inlet → expand hi until outlet(hi)>=Po
    lo=Po; hi=Po+10
    while hi < Po+2000:
        if calc_outlet(hi,Tc,L,D,Q,g) >= Po: break
        hi += 10
    # bisect: outlet(mid)<Po → Pi too low → lo=mid; else hi=mid
    for _ in range(60):
        mid=(lo+hi)/2; vm=calc_outlet(mid,Tc,L,D,Q,g)
        if abs(vm-Po)<0.005: return mid
        if vm < Po: lo=mid
        else:       hi=mid
    return (lo+hi)/2

def pipe_opts(D,P,gas):
    if gas=="O2": return [['1" (S14)','Required for O2']],'1" (S14)'
    rows,rec=[],None
    if D<=4.0:
        if P<=200:    rows=[['1/4" (S6)','up to 200 bar'],['1/4" (S9)','up to 1379 bar'],['1/4" (S12)','above 1379 bar']]; rec='1/4" (S6)'
        elif P<=1379: rows=[['1/4" (S9)','up to 1379 bar'],['1/4" (S12)','above 1379 bar']]; rec='1/4" (S9)'
        else:         rows=[['1/4" (S12)','above 1379 bar']]; rec='1/4" (S12)'
    elif D<=7.0:
        if P<=140: rows=[['3/8" (S16)','up to 140 bar'],['3/8" (S9)','above 140 bar']]; rec='3/8" (S16)'
        else:      rows=[['3/8" (S9)','above 140 bar']]; rec='3/8" (S9)'
    elif D<=21.0:
        if P<=20: rows=[['3/4" (S15)','up to 20 bar'],['1" (S14)','above 20 bar']]; rec='3/4" (S15)'
        else:     rows=[['1" (S14)','above 20 bar']]; rec='1" (S14)'
    else: rows=[['Special piping','Outside standard range']]; rec='Special piping'
    return rows,rec

# ── Session state ──────────────────────────────────────────────────────────
for k,v in dict(result=None,error=False,specs=[],
                gas_sel="N2",calc_sel="Pipe Diameter (mm)",field_vals={}).items():
    if k not in st.session_state: st.session_state[k]=v
if not st.session_state.field_vals:
    st.session_state.field_vals={k:str(v) for k,v in DEFAULTS.items()}

def fkey(f):
    return f.replace(" ","_").replace("(","").replace(")","").replace("/","").replace("°","deg")

# ── Query-param callback ───────────────────────────────────────────────────
qp=st.query_params
if qp.get("action")=="calc":
    try:
        gas=qp.get("gas","N2"); calc=qp.get("calc","Pipe Diameter (mm)")
        def gv(f): return float(qp.get(fkey(f), DEFAULTS.get(f,0.0)))
        T=gv("Temperature (°C)"); Pi=gv("Inlet Pressure (bar)")
        Po=gv("Outlet Pressure (bar)"); L=gv("Pipe Length (m)")
        D=gv("Pipe Diameter (mm)"); Q=gv("Flow Rate (LPM)")

        if   calc=="Pipe Diameter (mm)":
            r=calc_diameter(Pi,Po,T,L,Q,gas); sp,rec=pipe_opts(r,max(Pi,Po),gas)
            msg=f"Required Diameter: {r:.2f} mm<br>Recommended: {rec}"
        elif calc=="Flow Rate (LPM)":
            r=calc_flow(Pi,Po,T,L,D,gas); sp,rec=pipe_opts(D,max(Pi,Po),gas)
            msg=f"Maximum Flow Rate: {r:.1f} L/min<br>Recommended: {rec}"
        elif calc=="Pipe Length (m)":
            r=calc_length(Pi,Po,T,D,Q,gas); sp,rec=pipe_opts(D,max(Pi,Po),gas)
            msg=f"Maximum Pipe Length: {r:.1f} m<br>Recommended: {rec}"
        elif calc=="Inlet Pressure (bar)":
            r=calc_inlet(Po,T,L,D,Q,gas); sp,rec=pipe_opts(D,max(r,Po),gas)
            msg=f"Required Inlet Pressure: {r:.2f} bar<br>Recommended: {rec}"
        elif calc=="Outlet Pressure (bar)":
            r=calc_outlet(Pi,T,L,D,Q,gas); sp,rec=pipe_opts(D,Pi,gas)
            msg=f"Estimated Outlet Pressure: {r:.2f} bar<br>Recommended: {rec}"
        else:
            msg="Unsupported calculation type."; sp=[]

        st.session_state.update(result=msg, error=False, specs=sp,
                                gas_sel=gas, calc_sel=calc)
        for f in CALC_FIELDS.get(calc,[]):
            if fkey(f) in qp: st.session_state.field_vals[f]=qp[fkey(f)]
    except Exception as e:
        st.session_state.update(result=f"Error: {e}", error=True, specs=[])
    st.query_params.clear()

# ── Load HUD ───────────────────────────────────────────────────────────────
def load_hud()->str:
    p=Path(__file__).resolve().parent/"hud_background.html"
    if not p.exists(): return ""
    return "data:text/html;base64,"+base64.b64encode(p.read_bytes()).decode()

hud_src=load_hud()

# ── Build HTML fragments ───────────────────────────────────────────────────
gas_opts_html = "\n".join(
    f'<option value="{k}"{"selected" if k==st.session_state.gas_sel else ""}>{GAS_LABEL[k]}</option>'
    for k in GAS_DATA)

calc_opts_html = "\n".join(
    f'<option value="{k}"{"selected" if k==st.session_state.calc_sel else ""}>{k}</option>'
    for k in CALC_FIELDS)

cur_fields = CALC_FIELDS[st.session_state.calc_sel]
left_fields  = cur_fields[:4]
right_fields = cur_fields[4:]

def render_field_html(f):
    v=st.session_state.field_vals.get(f, DEFAULTS.get(f,0))
    k=fkey(f)
    return f'<div class="frow"><label>{f}</label><input type="number" id="{k}" value="{v}" step="any"></div>'

left_html  = "".join(render_field_html(f) for f in left_fields)
right_html = "".join(render_field_html(f) for f in right_fields)

def make_result_html():
    if not st.session_state.result: return ""
    if st.session_state.error:
        return f'<div class="rbox err">⚠ {st.session_state.result}</div>'
    rows="".join(f"<tr><td>{i}</td><td>{r[0]}</td><td>{r[1]}</td></tr>"
                 for i,r in enumerate(st.session_state.specs))
    tbl=(f'<p class="sp-title">Possible pipe specifications:</p>'
         f'<table class="stbl"><thead><tr><th>#</th><th>Pipe Spec</th><th>Details</th></tr></thead>'
         f'<tbody>{rows}</tbody></table>') if st.session_state.specs else ""
    return f'<div class="rbox ok">{st.session_state.result}</div>{tbl}'

result_html=make_result_html()

# ── JS field/calc definitions ──────────────────────────────────────────────
CF_JS="""{
  "Pipe Diameter (mm)":    ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  "Flow Rate (LPM)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  "Pipe Length (m)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Inlet Pressure (bar)":  ["Temperature (°C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Outlet Pressure (bar)": ["Temperature (°C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"]
}"""
DEF_JS="""{
  "Temperature (°C)":25,"Inlet Pressure (bar)":100,
  "Outlet Pressure (bar)":10,"Pipe Length (m)":10,
  "Pipe Diameter (mm)":10,"Flow Rate (LPM)":100
}"""

# ── Full page HTML ─────────────────────────────────────────────────────────
page_html=f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
/* ══ Reset ══ */
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
html,body{{
  width:100%;height:100%;
  overflow:hidden;          /* page never scrolls */
  background:#000d0d;
  font-family:'Courier New',monospace;
  color:#00e5cc;
}}

/* ══ HUD iframe: true full-screen behind everything ══ */
#hud{{
  position:fixed;inset:0;
  width:100%;height:100%;
  border:none;pointer-events:none;
  z-index:1;
}}

/* ══════════════════════════════════════════════════
   LAYOUT
   Columns: 2vw | 34vw | 1fr | 34vw | 2vw
   Rows:    17vh (HUD top) | auto (calc) | 30vh (HUD bottom)
   ══════════════════════════════════════════════════ */
#grid{{
  position:fixed;inset:0;
  display:grid;
  grid-template-columns: 2vw 34vw 1fr 34vw 2vw;
  grid-template-rows:    17vh auto 30vh;
  pointer-events:none;
  z-index:500;
}}

/* Title — row 1, cols 2-4 */
#title-cell{{
  grid-row:1; grid-column:2/5;
  display:flex;
  align-items:flex-end;
  justify-content:center;
  padding-bottom:8px;
  pointer-events:none;
}}
#title-inner{{
  text-align:center;
  line-height:1.35;
}}
#title-inner h1{{
  color:#00ffee;
  text-shadow:0 0 16px #00e5cc,0 0 32px rgba(0,229,204,0.35);
  font-size:clamp(14px,1.6vw,22px);
  letter-spacing:5px;
  text-transform:uppercase;
  font-family:'Courier New',monospace;
}}
#title-inner h1 .tri{{
  color:#00e5cc;
  margin-right:10px;
}}

/* Left panel — row 2, col 2 */
#left{{
  grid-row:2; grid-column:2;
  pointer-events:all;
  background:#000000;
  border:1px solid rgba(0,229,204,0.35);
  border-radius:3px;
  padding:14px 16px 16px;
  overflow-y:auto;
  overflow-x:hidden;
  /* height = fit content; scrolls if taller than available space */
  align-self:start;
  max-height:calc(53vh);    /* cap so it never overflows into HUD bottom */
  scrollbar-width:thin;
  scrollbar-color:#00e5cc #001813;
}}
#left::-webkit-scrollbar{{width:4px}}
#left::-webkit-scrollbar-thumb{{background:#00e5cc;border-radius:2px}}

/* Right panel — row 2, col 4 */
#right{{
  grid-row:2; grid-column:4;
  pointer-events:all;
  background:#000000;
  border:1px solid rgba(0,229,204,0.35);
  border-radius:3px;
  padding:14px 16px 16px;
  overflow-y:auto;
  overflow-x:hidden;
  align-self:start;
  max-height:calc(53vh);
  scrollbar-width:thin;
  scrollbar-color:#00e5cc #001813;
}}
#right::-webkit-scrollbar{{width:4px}}
#right::-webkit-scrollbar-thumb{{background:#00e5cc;border-radius:2px}}

/* ══ Shared panel content styles ══ */

/* select rows */
.sr{{display:flex;flex-direction:column;gap:3px;margin-bottom:11px}}
.sr>label{{
  font-size:9.5px;letter-spacing:1.5px;text-transform:uppercase;
  color:rgba(0,229,204,0.6);
}}
.sr select{{
  background:#000;color:#00ffee;
  border:1px solid rgba(0,229,204,0.4);border-radius:2px;
  padding:6px 26px 6px 8px;
  font-family:'Courier New',monospace;font-size:11.5px;
  width:100%;cursor:pointer;outline:none;
  -webkit-appearance:none;appearance:none;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%2300e5cc'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 9px center;
}}
.sr select:focus{{border-color:#00e5cc;box-shadow:0 0 7px rgba(0,229,204,0.35)}}
.sr select option{{background:#000d0d;color:#00e5cc}}

/* divider */
.dv{{border:none;border-top:1px solid rgba(0,229,204,0.18);margin:10px 0 12px}}

/* field row */
.frow{{margin-bottom:9px}}
.frow>label{{
  display:block;font-size:10px;
  color:rgba(0,229,204,0.75);
  margin-bottom:4px;letter-spacing:0.4px;
}}
.frow input{{
  width:100%;
  background:#000;color:#00ffee;
  border:1px solid rgba(0,229,204,0.4);border-radius:2px;
  padding:6px 8px;
  font-family:'Courier New',monospace;font-size:12px;
  outline:none;
}}
.frow input:focus{{border-color:#00e5cc;box-shadow:0 0 6px rgba(0,229,204,0.3)}}
.frow input::-webkit-inner-spin-button,
.frow input::-webkit-outer-spin-button{{-webkit-appearance:none;margin:0}}
.frow input[type=number]{{-moz-appearance:textfield}}

/* caption */
.cap{{
  font-size:9px;color:rgba(0,229,204,0.4);
  letter-spacing:1px;margin:4px 0 13px;
}}

/* CALCULATE button */
#cbtn{{
  display:block;width:100%;
  background:rgba(0,229,204,0.07);color:#00ffee;
  border:1.5px solid #00e5cc;border-radius:3px;
  padding:9px 0;
  font-family:'Courier New',monospace;
  font-size:12px;letter-spacing:4px;
  cursor:pointer;text-transform:uppercase;
  transition:background 0.2s,box-shadow 0.2s;
  margin-bottom:14px;
}}
#cbtn:hover{{background:rgba(0,229,204,0.2);box-shadow:0 0 18px rgba(0,229,204,0.5)}}
#cbtn:active{{background:rgba(0,229,204,0.32)}}

/* result box — green tint as in screenshot */
.rbox{{
  padding:10px 13px;border-radius:3px;
  font-size:11.5px;line-height:1.9;word-break:break-word;
  margin-bottom:12px;
}}
.rbox.ok{{
  background:rgba(0,80,40,0.85);
  border:1px solid rgba(0,229,100,0.5);
  border-left:3px solid #00e87a;
  color:#00ffcc;
  font-weight:bold;
}}
.rbox.err{{
  background:rgba(60,0,0,0.9);
  border:1px solid rgba(255,60,60,0.4);
  border-left:3px solid #ff4444;
  color:#ff9090;
}}

/* spec table */
.sp-title{{
  font-size:10px;letter-spacing:1px;
  color:#00e5cc;text-transform:uppercase;
  margin-bottom:6px;
}}
.stbl{{width:100%;border-collapse:collapse;font-size:10.5px}}
.stbl th{{
  color:#00e5cc;background:rgba(0,229,204,0.08);
  border-bottom:1px solid rgba(0,229,204,0.3);
  padding:5px 8px;text-align:left;
  letter-spacing:1px;font-size:9.5px;text-transform:uppercase;
}}
.stbl td{{
  color:rgba(0,229,204,0.85);
  border-bottom:1px solid rgba(0,229,204,0.1);
  padding:5px 8px;
}}
.stbl tr:last-child td{{border-bottom:none}}
</style>
</head>
<body>

<!-- HUD: full-screen animated background (z-index 1) -->
{"<iframe id='hud' src='" + hud_src + "' sandbox='allow-scripts'></iframe>" if hud_src else ""}

<!-- Grid layout (z-index 500, pointer-events:none on grid; all on panels) -->
<div id="grid">

  <!-- TITLE (row1, cols2-4) -->
  <div id="title-cell">
    <div id="title-inner">
      <h1><span class="tri">⟁</span>HIGH-PRESSURE GAS FLOW<br>CALCULATOR</h1>
    </div>
  </div>

  <!-- LEFT PANEL (row2, col2) -->
  <div id="left">
    <div class="sr">
      <label>Gas type:</label>
      <select id="gasSelect" onchange="rebuildAll()">{gas_opts_html}</select>
    </div>
    <div class="sr">
      <label>Calculation type:</label>
      <select id="calcSelect" onchange="rebuildAll()">{calc_opts_html}</select>
    </div>
    <hr class="dv">
    <div id="left-fields">{left_html}</div>
  </div>

  <!-- RIGHT PANEL (row2, col4) -->
  <div id="right">
    <div id="right-fields">{right_html}</div>
    <div class="cap">Friction factor (f) = 0.02 · fixed constant</div>
    <button id="cbtn" onclick="doCalc()">▶ &nbsp;Calculate</button>
    <div id="result-area">{result_html}</div>
  </div>

</div><!-- #grid -->

<script>
const CF={CF_JS};
const DEF={DEF_JS};

function toKey(f){{
  return f.replace(/ /g,'_').replace(/[()]/g,'').replace(/\\//g,'').replace(/°/g,'deg');
}}

function renderField(f){{
  const k=toKey(f), v=DEF[f]??0;
  return `<div class="frow"><label>${{f}}</label><input type="number" id="${{k}}" value="${{v}}" step="any"></div>`;
}}

function rebuildAll(){{
  const calc=document.getElementById('calcSelect').value;
  const all=CF[calc]||[];
  document.getElementById('left-fields').innerHTML =all.slice(0,4).map(renderField).join('');
  document.getElementById('right-fields').innerHTML=all.slice(4).map(renderField).join('');
  document.getElementById('result-area').innerHTML ='';
}}

function doCalc(){{
  const gas =document.getElementById('gasSelect').value;
  const calc=document.getElementById('calcSelect').value;
  const p   =new URLSearchParams();
  p.set('action','calc'); p.set('gas',gas); p.set('calc',calc);
  (CF[calc]||[]).forEach(f=>{{
    const el=document.getElementById(toKey(f));
    p.set(toKey(f), el?el.value:(DEF[f]??0));
  }});
  window.parent.location.search='?'+p.toString();
}}
</script>
</body>
</html>"""

# ── Streamlit: hide chrome, expand iframe to full viewport ─────────────────
st.markdown("""
<style>
  #MainMenu,footer,header,
  [data-testid="stToolbar"],
  [data-testid="stDecoration"],
  [data-testid="stStatusWidget"],
  [data-testid="collapsedControl"]{{display:none!important}}
  .stApp,[data-testid="stAppViewContainer"],section.main{{
    background:#000d0d!important;padding:0!important;margin:0!important}}
  section.main .block-container{{
    padding:0!important;margin:0!important;max-width:100vw!important}}
  iframe[title="streamlit_components_v1.html_v1"]{{
    position:fixed!important;inset:0!important;
    width:100vw!important;height:100vh!important;
    border:none!important;z-index:100!important}}
</style>
""", unsafe_allow_html=True)

components.html(page_html, height=800, scrolling=False)
