"""
High-Pressure Gas Flow Calculator
Layout matches screenshot:
  • Full-screen HUD animated background
  • Title centered at top-center (over dark backdrop, within HUD space)
  • Left panel  : Gas type | Calc type | first 4 input fields
  • Right panel : last input field | Calculate btn | result | spec table
  • HUD graphics visible: top, bottom, far-left, far-right strips
  • Neither panel overlaps the HUD border zones

Run: streamlit run app.py
(hud_background.html must be in the same folder)
"""

import math
import base64
import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(page_title="Gas Flow Calculator", page_icon="💨",
                   layout="wide", initial_sidebar_state="collapsed")

# ─────────────────────────────────────────────────────────────────────────────
# DATA
# ─────────────────────────────────────────────────────────────────────────────
GAS_DATA = {
    "N2":0.028013,"O2":0.031999,"Ar":0.039948,"CO2":0.04401,
    "He":0.0040026,"H2":0.002016,"CH4":0.01604,"C2H2":0.02604,
    "Forming Gas 1":0.03881,"Forming Gas 2":0.02671,"Air":0.02897,
}
GAS_LABEL = {
    "N2":"N2 (Nitrogen)","O2":"O2 (Oxygen)","Ar":"Ar (Argon)",
    "CO2":"CO2 (Carbon Dioxide)","He":"He (Helium)","H2":"H2 (Hydrogen)",
    "CH4":"CH4 (Methane)","C2H2":"C2H2 (Acetylene)",
    "Forming Gas 1":"H2 3% + Ar 97% (FG1)","Forming Gas 2":"H2 5% + N2 95% (FG2)",
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

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICS
# ─────────────────────────────────────────────────────────────────────────────
def _rho(P,T,g):      return (P*GAS_DATA[g])/(8.314*T)
def _ra(Pi,Po,T,g):   return (_rho(Pi*1e5,T,g)+_rho(Po*1e5,T,g))/2
F=0.02

def calc_diameter(Pi,Po,Tc,L,Q,g):
    T=Tc+273.15;Qs=Q/6e4;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet pressure must exceed outlet pressure.")
    return ((F*L*8*_ra(Pi,Po,T,g)*Qs**2)/(math.pi**2*dP))**0.2*1000

def calc_flow(Pi,Po,Tc,L,D,g):
    T=Tc+273.15;Dm=D/1000;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet pressure must exceed outlet pressure.")
    return math.sqrt((dP*math.pi**2*Dm**5)/(8*F*L*_ra(Pi,Po,T,g)))*6e4

def calc_length(Pi,Po,Tc,D,Q,g):
    T=Tc+273.15;Dm=D/1000;Qs=Q/6e4;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet pressure must exceed outlet pressure.")
    return (dP*math.pi**2*Dm**5)/(8*F*_ra(Pi,Po,T,g)*Qs**2)

def calc_outlet(Pi,Tc,L,D,Q,g):
    T=Tc+273.15;Dm=D/1000;Qs=Q/6e4;Pg=Pi
    for _ in range(20):
        ra=(_rho(Pi*1e5,T,g)+_rho(Pg*1e5,T,g))/2
        dP=(8*F*L*ra*Qs**2)/(math.pi**2*Dm**5)
        Pn=Pi-dP/1e5
        if abs(Pn-Pg)<0.001: return max(Pn,0.0)
        Pg=Pn
    return max(Pg,0.0)

def calc_inlet(Po,Tc,L,D,Q,g):
    lo,hi=Po,Po+100
    while hi<Po+2000:
        if calc_outlet(hi,Tc,L,D,Q,g)<=Po: break
        hi+=100
    for _ in range(50):
        mid=(lo+hi)/2;vm=calc_outlet(mid,Tc,L,D,Q,g)
        if abs(vm-Po)<0.01: return mid
        lo,hi=(mid,hi) if vm>Po else (lo,mid)
    return (lo+hi)/2

def pipe_opts(D,P,gas):
    if gas=="O2": return [['1" (S14)','Required for O2']],'1" (S14)'
    rows,rec=[],None
    if D<=4.0:
        if P<=200:    rows=[['1/4" (S6)','up to 200 bar'],['1/4" (S9)','up to 1379 bar'],['1/4" (S12)','above 1379 bar']];rec='1/4" (S6)'
        elif P<=1379: rows=[['1/4" (S9)','up to 1379 bar'],['1/4" (S12)','above 1379 bar']];rec='1/4" (S9)'
        else:         rows=[['1/4" (S12)','above 1379 bar']];rec='1/4" (S12)'
    elif D<=7.0:
        if P<=140: rows=[['3/8" (S16)','up to 140 bar'],['3/8" (S9)','above 140 bar']];rec='3/8" (S16)'
        else:      rows=[['3/8" (S9)','above 140 bar']];rec='3/8" (S9)'
    elif D<=21.0:
        if P<=20: rows=[['3/4" (S15)','up to 20 bar'],['1" (S14)','above 20 bar']];rec='3/4" (S15)'
        else:     rows=[['1" (S14)','above 20 bar']];rec='1" (S14)'
    else: rows=[['Special piping','Outside standard range']];rec='Special piping'
    return rows,rec

# ─────────────────────────────────────────────────────────────────────────────
# SESSION STATE
# ─────────────────────────────────────────────────────────────────────────────
for k,v in dict(result=None,error=False,specs=[],
                gas_sel="N2",calc_sel="Pipe Diameter (mm)",field_vals={}).items():
    if k not in st.session_state: st.session_state[k]=v
if not st.session_state.field_vals:
    st.session_state.field_vals={k:str(v) for k,v in DEFAULTS.items()}

def fkey(f):
    return f.replace(" ","_").replace("(","").replace(")","").replace("/","").replace("°","deg")

# ─────────────────────────────────────────────────────────────────────────────
# QUERY PARAM CALLBACK
# ─────────────────────────────────────────────────────────────────────────────
qp=st.query_params
if qp.get("action")=="calc":
    try:
        gas=qp.get("gas","N2"); calc=qp.get("calc","Pipe Diameter (mm)")
        def gv(f): return float(qp.get(fkey(f),DEFAULTS.get(f,0.0)))
        T,Pi,Po,L,D,Q = gv("Temperature (°C)"),gv("Inlet Pressure (bar)"),gv("Outlet Pressure (bar)"),gv("Pipe Length (m)"),gv("Pipe Diameter (mm)"),gv("Flow Rate (LPM)")
        if   calc=="Pipe Diameter (mm)":    r=calc_diameter(Pi,Po,T,L,Q,gas);sp,rec=pipe_opts(r,max(Pi,Po),gas);msg=f"Required Diameter: {r:.2f} mm<br>Recommended: {rec}"
        elif calc=="Flow Rate (LPM)":       r=calc_flow(Pi,Po,T,L,D,gas);sp,rec=pipe_opts(D,max(Pi,Po),gas);msg=f"Maximum Flow Rate: {r:.1f} L/min<br>Recommended: {rec}"
        elif calc=="Pipe Length (m)":       r=calc_length(Pi,Po,T,D,Q,gas);sp,rec=pipe_opts(D,max(Pi,Po),gas);msg=f"Maximum Pipe Length: {r:.1f} m<br>Recommended: {rec}"
        elif calc=="Inlet Pressure (bar)":  r=calc_inlet(Po,T,L,D,Q,gas);sp,rec=pipe_opts(D,max(r,Po),gas);msg=f"Required Inlet Pressure: {r:.2f} bar<br>Recommended: {rec}"
        elif calc=="Outlet Pressure (bar)": r=calc_outlet(Pi,T,L,D,Q,gas);sp,rec=pipe_opts(D,Pi,gas);msg=f"Estimated Outlet Pressure: {r:.2f} bar<br>Recommended: {rec}"
        else: msg="Unsupported type.";sp=[]
        st.session_state.update(result=msg,error=False,specs=sp,gas_sel=gas,calc_sel=calc)
        for f in CALC_FIELDS.get(calc,[]):
            if fkey(f) in qp: st.session_state.field_vals[f]=qp[fkey(f)]
    except Exception as e:
        st.session_state.update(result=f"Error: {e}",error=True,specs=[])
    st.query_params.clear()

# ─────────────────────────────────────────────────────────────────────────────
# LOAD HUD
# ─────────────────────────────────────────────────────────────────────────────
def load_hud()->str:
    p=Path(__file__).resolve().parent/"hud_background.html"
    if not p.exists(): return ""
    return "data:text/html;base64,"+base64.b64encode(p.read_bytes()).decode()

hud_src=load_hud()

# ─────────────────────────────────────────────────────────────────────────────
# BUILD HTML FRAGMENTS
# ─────────────────────────────────────────────────────────────────────────────
gas_opts="\n".join(
    f'<option value="{k}"{"selected" if k==st.session_state.gas_sel else ""}>{GAS_LABEL[k]}</option>'
    for k in GAS_DATA)

calc_opts="\n".join(
    f'<option value="{k}"{"selected" if k==st.session_state.calc_sel else ""}>{k}</option>'
    for k in CALC_FIELDS)

# Split fields: first 4 go left, last 1 goes right
cur_fields = CALC_FIELDS[st.session_state.calc_sel]
left_fields  = cur_fields[:4]
right_fields = cur_fields[4:]

def render_field(f):
    v = st.session_state.field_vals.get(f, DEFAULTS.get(f, 0))
    k = fkey(f)
    return f"""<div class="frow">
      <label>{f}</label>
      <div class="input-wrap">
        <input type="number" id="{k}" value="{v}" step="any">
        <button class="spin-btn" onclick="adj('{k}',-1)">−</button>
        <button class="spin-btn" onclick="adj('{k}',1)">+</button>
      </div>
    </div>"""

left_fields_html  = "".join(render_field(f) for f in left_fields)
right_fields_html = "".join(render_field(f) for f in right_fields)

# Result
def make_result():
    if not st.session_state.result: return ""
    if st.session_state.error:
        return f'<div class="rbox err">⚠ {st.session_state.result}</div>'
    spec_rows="".join(f"<tr><td>{i}</td><td>{r[0]}</td><td>{r[1]}</td></tr>"
                      for i,r in enumerate(st.session_state.specs))
    spec_tbl=(f'<p class="spec-title">Possible pipe specifications:</p>'
              f'<table class="stbl"><thead><tr><th></th><th>Pipe Spec</th><th>Details</th></tr></thead>'
              f'<tbody>{spec_rows}</tbody></table>') if st.session_state.specs else ""
    return f'<div class="rbox ok">{st.session_state.result}</div>{spec_tbl}'

result_html=make_result()

# ─────────────────────────────────────────────────────────────────────────────
# JAVASCRIPT field definitions (for dynamic rebuild)
# ─────────────────────────────────────────────────────────────────────────────
cf_js = """{
  "Pipe Diameter (mm)":    ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  "Flow Rate (LPM)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  "Pipe Length (m)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Inlet Pressure (bar)":  ["Temperature (°C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Outlet Pressure (bar)": ["Temperature (°C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"]
}"""

def_js = """{
  "Temperature (°C)":25,"Inlet Pressure (bar)":100,
  "Outlet Pressure (bar)":10,"Pipe Length (m)":10,
  "Pipe Diameter (mm)":10,"Flow Rate (LPM)":100
}"""

# ─────────────────────────────────────────────────────────────────────────────
# FULL PAGE HTML
# ─────────────────────────────────────────────────────────────────────────────
page_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
html,body{{width:100%;height:100%;overflow:hidden;background:#000d0d;font-family:'Courier New',monospace;color:#00e5cc}}

/* ══ HUD full-screen behind everything ══ */
#hud{{position:fixed;inset:0;width:100%;height:100%;border:none;pointer-events:none;z-index:1}}

/* ══ Layout grid ══
   Rows:
     row1: 17vh  — HUD top strip (top panels, waveforms)
     row2: auto  — calculator area (title + left/right panels)
     row3: 30vh  — HUD bottom strip (binary, access, tracking, analyzing)
   Cols:
     col1: 2vw   — left HUD margin
     col2: 35vw  — left calc panel
     col3: 1fr   — gap / center HUD
     col4: 35vw  — right calc panel
     col5: 2vw   — right HUD margin
*/
#layout{{
  position:fixed;inset:0;
  display:grid;
  grid-template-rows:16vh 1fr 32vh;
  grid-template-columns:4vw 32vw 1fr 32vw 4vw;
  z-index:500;
  pointer-events:none;
}}

/* ══ Title bar: spans cols 2-4, row 1 bottom ══ */
#title-bar{{
  grid-row:1;grid-column:2/5;
  display:flex;align-items:flex-end;justify-content:center;
  padding-bottom:10px;
  pointer-events:none;
}}
#title-inner{{
  background:rgba(0,0,0,0.0);
  padding:6px 30px 4px;
  text-align:center;
}}
#title-inner h1{{
  color:#00ffee;
  text-shadow:0 0 16px #00e5cc,0 0 30px rgba(0,229,204,0.4);
  font-size:22px;
  letter-spacing:6px;
  text-transform:uppercase;
  font-family:'Courier New',monospace;
  line-height:1.3;
}}
#title-inner h1 span{{color:#00e5cc;margin-right:8px}}

/* ══ Left panel: row 2, col 2 ══ */
#left-panel{{
  grid-row:2;grid-column:2;
  pointer-events:all;
  background:rgba(0,0,0,0.85);
  border:1px solid rgba(0,229,204,0.35);
  border-radius:3px;
  padding:16px 18px;
  overflow-y:auto;
  overflow-x:hidden;
  scrollbar-width:thin;scrollbar-color:#00e5cc #001a18;
  align-self:stretch;
  margin:4px 0;
}}
#left-panel::-webkit-scrollbar{{width:4px}}
#left-panel::-webkit-scrollbar-thumb{{background:#00e5cc;border-radius:2px}}

/* ══ Right panel: row 2, col 4 ══ */
#right-panel{{
  grid-row:2;grid-column:4;
  pointer-events:all;
  background:rgba(0,0,0,0.85);
  border:1px solid rgba(0,229,204,0.35);
  border-radius:3px;
  padding:16px 18px;
  overflow-y:auto;
  overflow-x:hidden;
  scrollbar-width:thin;scrollbar-color:#00e5cc #001a18;
  align-self:stretch;
  margin:4px 0;
}}
#right-panel::-webkit-scrollbar{{width:4px}}
#right-panel::-webkit-scrollbar-thumb{{background:#00e5cc;border-radius:2px}}

/* ══ Shared panel styles ══ */
.panel-section-label{{
  font-size:9px;letter-spacing:2px;text-transform:uppercase;
  color:rgba(0,229,204,0.5);margin-bottom:6px;
}}

.srow{{display:flex;flex-direction:column;gap:3px;margin-bottom:12px}}
.srow label{{font-size:10px;letter-spacing:1px;color:rgba(0,229,204,0.7);text-transform:uppercase}}
.srow select{{
  background:#000;color:#00ffee;
  border:1px solid rgba(0,229,204,0.4);border-radius:2px;
  padding:6px 28px 6px 8px;
  font-family:'Courier New',monospace;font-size:11px;
  width:100%;cursor:pointer;outline:none;
  appearance:none;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%2300e5cc'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 10px center;
}}
.srow select:focus{{border-color:#00e5cc}}
.srow select option{{background:#000d0d;color:#00e5cc}}

.hdiv{{border:none;border-top:1px solid rgba(0,229,204,0.15);margin:12px 0}}

/* input field row */
.frow{{margin-bottom:10px}}
.frow label{{display:block;font-size:10px;color:rgba(0,229,204,0.75);margin-bottom:4px;letter-spacing:0.5px}}
.input-wrap{{display:flex;align-items:stretch;gap:0}}
.input-wrap input{{
  flex:1;min-width:0;
  background:#000;color:#00ffee;
  border:1px solid rgba(0,229,204,0.4);
  border-right:none;
  border-radius:2px 0 0 2px;
  padding:6px 8px;
  font-family:'Courier New',monospace;font-size:12px;outline:none;
}}
.input-wrap input:focus{{border-color:#00e5cc;box-shadow:0 0 6px rgba(0,229,204,0.3)}}
.input-wrap input::-webkit-inner-spin-button,
.input-wrap input::-webkit-outer-spin-button{{-webkit-appearance:none;margin:0}}
.input-wrap input[type=number]{{-moz-appearance:textfield}}
.spin-btn{{
  background:rgba(0,229,204,0.1);color:#00e5cc;
  border:1px solid rgba(0,229,204,0.4);
  padding:0 10px;cursor:pointer;
  font-size:14px;font-family:'Courier New',monospace;
  transition:background 0.15s;
  line-height:1;
}}
.spin-btn:last-child{{border-radius:0 2px 2px 0}}
.spin-btn:first-of-type{{border-left:none}}
.spin-btn:hover{{background:rgba(0,229,204,0.25)}}

.caption{{font-size:9px;color:rgba(0,229,204,0.4);letter-spacing:1px;margin:4px 0 14px}}

/* calculate button */
#calc-btn{{
  display:block;width:100%;
  background:rgba(0,229,204,0.08);color:#00ffee;
  border:1.5px solid #00e5cc;border-radius:3px;
  padding:9px 0;
  font-family:'Courier New',monospace;font-size:12px;letter-spacing:4px;
  cursor:pointer;text-transform:uppercase;
  transition:background 0.2s,box-shadow 0.2s;
  margin-bottom:14px;
}}
#calc-btn:hover{{background:rgba(0,229,204,0.2);box-shadow:0 0 18px rgba(0,229,204,0.5)}}
#calc-btn:active{{background:rgba(0,229,204,0.32)}}

/* result */
.rbox{{
  padding:10px 12px;border-radius:3px;
  font-size:11.5px;line-height:1.9;word-break:break-word;
  margin-bottom:12px;
}}
.rbox.ok {{background:rgba(0,40,30,0.95);border:1px solid rgba(0,229,204,0.3);border-left:3px solid #00e5cc;color:#00ffee}}
.rbox.err{{background:rgba(40,0,0,0.95);border:1px solid rgba(255,60,60,0.3);border-left:3px solid #ff4444;color:#ff9090}}

.spec-title{{font-size:10px;letter-spacing:1px;color:#00e5cc;margin-bottom:6px;text-transform:uppercase}}

.stbl{{width:100%;border-collapse:collapse;font-size:10.5px}}
.stbl th{{
  color:#00e5cc;background:rgba(0,229,204,0.08);
  border-bottom:1px solid rgba(0,229,204,0.3);
  padding:5px 8px;text-align:left;letter-spacing:1px;font-size:9.5px;text-transform:uppercase;
}}
.stbl td{{color:rgba(0,229,204,0.82);border-bottom:1px solid rgba(0,229,204,0.1);padding:5px 8px}}
.stbl tr:last-child td{{border-bottom:none}}
</style>
</head>
<body>

<!-- HUD: full-screen animated background -->
{"<iframe id='hud' src='" + hud_src + "' sandbox='allow-scripts'></iframe>" if hud_src else ""}

<!-- Main layout grid -->
<div id="layout">

  <!-- Title (row1, cols 2-4) -->
  <div id="title-bar">
    <div id="title-inner">
      <h1><span>⟁</span>HIGH-PRESSURE GAS FLOW<br>CALCULATOR</h1>
    </div>
  </div>

  <!-- LEFT PANEL (row2, col2) -->
  <div id="left-panel">
    <div class="srow">
      <label>Gas type:</label>
      <select id="gasSelect" onchange="rebuildAll()">{gas_opts}</select>
    </div>

    <div class="srow">
      <label>Calculation type:</label>
      <select id="calcSelect" onchange="rebuildAll()">{calc_opts}</select>
    </div>

    <hr class="hdiv">

    <div id="left-fields">{left_fields_html}</div>
  </div>

  <!-- RIGHT PANEL (row2, col4) -->
  <div id="right-panel">
    <div id="right-fields">{right_fields_html}</div>

    <div class="caption">Friction factor (f) = 0.02 · fixed constant</div>

    <button id="calc-btn" onclick="doCalc()">▶ &nbsp;Calculate</button>

    <div id="result-area">{result_html}</div>
  </div>

</div><!-- #layout -->

<script>
const CF = {cf_js};
const DEF = {def_js};

function toKey(f){{
  return f.replace(/ /g,'_').replace(/[()]/g,'').replace(/\//g,'').replace(/°/g,'deg');
}}

function renderField(f){{
  const k=toKey(f), v=DEF[f]??0;
  return `<div class="frow">
    <label>${{f}}</label>
    <div class="input-wrap">
      <input type="number" id="${{k}}" value="${{v}}" step="any">
      <button class="spin-btn" onclick="adj('${{k}}',-1)">−</button>
      <button class="spin-btn" onclick="adj('${{k}}',1)">+</button>
    </div>
  </div>`;
}}

function adj(id, delta){{
  const el=document.getElementById(id);
  if(!el) return;
  const v=parseFloat(el.value)||0;
  // step = 1 for most, 0.1 for temperature
  const step = id.includes('Temperature') ? 0.5 : 1;
  el.value = Math.round((v + delta*step)*1000)/1000;
}}

function rebuildAll(){{
  const calc=document.getElementById('calcSelect').value;
  const fields=CF[calc]||[];
  const leftF=fields.slice(0,4);
  const rightF=fields.slice(4);
  document.getElementById('left-fields').innerHTML  = leftF.map(renderField).join('');
  document.getElementById('right-fields').innerHTML = rightF.map(renderField).join('');
  document.getElementById('result-area').innerHTML  = '';
}}

function doCalc(){{
  const gas =document.getElementById('gasSelect').value;
  const calc=document.getElementById('calcSelect').value;
  const p   =new URLSearchParams();
  p.set('action','calc'); p.set('gas',gas); p.set('calc',calc);
  (CF[calc]||[]).forEach(f=>{{
    const el=document.getElementById(toKey(f));
    p.set(toKey(f), el ? el.value : (DEF[f]??0));
  }});
  window.parent.location.search='?'+p.toString();
}}
</script>
</body>
</html>"""

# ─────────────────────────────────────────────────────────────────────────────
# STREAMLIT OUTPUT
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
  #MainMenu,footer,header,
  [data-testid="stToolbar"],
  [data-testid="stDecoration"],
  [data-testid="stStatusWidget"],
  [data-testid="collapsedControl"]{{display:none!important}}
  .stApp,[data-testid="stAppViewContainer"],section.main{{
    background:#000d0d!important;padding:0!important;margin:0!important}}
  section.main .block-container{{padding:0!important;margin:0!important;max-width:100vw!important}}
  iframe[title="streamlit_components_v1.html_v1"]{{
    position:fixed!important;inset:0!important;
    width:100vw!important;height:100vh!important;
    border:none!important;z-index:100!important}}
</style>
""", unsafe_allow_html=True)

components.html(page_html, height=800, scrolling=False)
