import math
import re
import base64
import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
#  Page config
# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="High-Pressure Gas Flow Calculator",
    page_icon="💨",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ─────────────────────────────────────────────────────────────────────────────
#  Load HUD HTML → base64 data-URI so we can embed it as a true full-page iframe
# ─────────────────────────────────────────────────────────────────────────────
def load_hud_b64() -> str:
    p = Path(__file__).resolve().parent / "hud_background.html"
    if not p.exists():
        return ""
    raw = p.read_bytes()
    return "data:text/html;base64," + base64.b64encode(raw).decode()

hud_src = load_hud_b64()

# ─────────────────────────────────────────────────────────────────────────────
#  Domain data
# ─────────────────────────────────────────────────────────────────────────────
GAS_DATA = {
    "N2":0.028013,"O2":0.031999,"Ar":0.039948,"CO2":0.04401,
    "He":0.0040026,"H2":0.002016,"CH4":0.01604,"C2H2":0.02604,
    "Forming Gas 1":0.03881,"Forming Gas 2":0.02671,"Air":0.02897,
}
GAS_LABEL = {
    "N2":"N2 – Nitrogen","O2":"O2 – Oxygen","Ar":"Ar – Argon",
    "CO2":"CO2 – Carbon Dioxide","He":"He – Helium","H2":"H2 – Hydrogen",
    "CH4":"CH4 – Methane","C2H2":"C2H2 – Acetylene",
    "Forming Gas 1":"H2 3% + Ar 97%  (FG1)",
    "Forming Gas 2":"H2 5% + N2 95%  (FG2)",
    "Air":"Air – Dry Air",
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
#  Physics
# ─────────────────────────────────────────────────────────────────────────────
def _rho(P,T,g):  return (P*GAS_DATA[g])/(8.314*T)
def _ra(Pi,Po,T,g): return (_rho(Pi*1e5,T,g)+_rho(Po*1e5,T,g))/2

def calc_diameter(Pi,Po,Tc,L,Q,g):
    T=Tc+273.15;Qs=Q/6e4;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet must exceed outlet pressure.")
    return ((0.02*L*8*_ra(Pi,Po,T,g)*Qs**2)/(math.pi**2*dP))**0.2*1000

def calc_flow(Pi,Po,Tc,L,D,g):
    T=Tc+273.15;Dm=D/1000;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet must exceed outlet pressure.")
    return math.sqrt((dP*math.pi**2*Dm**5)/(8*0.02*L*_ra(Pi,Po,T,g)))*6e4

def calc_length(Pi,Po,Tc,D,Q,g):
    T=Tc+273.15;Dm=D/1000;Qs=Q/6e4;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet must exceed outlet pressure.")
    return (dP*math.pi**2*Dm**5)/(8*0.02*_ra(Pi,Po,T,g)*Qs**2)

def calc_outlet(Pi,Tc,L,D,Q,g):
    T=Tc+273.15;Dm=D/1000;Qs=Q/6e4;Pg=Pi
    for _ in range(20):
        ra=(_rho(Pi*1e5,T,g)+_rho(Pg*1e5,T,g))/2
        dP=(8*0.02*L*ra*Qs**2)/(math.pi**2*Dm**5)
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
        if vm>Po: lo=mid
        else: hi=mid
    return (lo+hi)/2

def pipe_opts(D,P,gas):
    if gas=="O2": return [['1" (S14)','Required for O2']],'1" (S14)'
    rows,rec=[],None
    if D<=4.0:
        if P<=200:    rows=[['1/4" (S6)','≤200 bar'],['1/4" (S9)','≤1379 bar'],['1/4" (S12)','>1379 bar']];rec='1/4" (S6)'
        elif P<=1379: rows=[['1/4" (S9)','≤1379 bar'],['1/4" (S12)','>1379 bar']];rec='1/4" (S9)'
        else:         rows=[['1/4" (S12)','>1379 bar']];rec='1/4" (S12)'
    elif D<=7.0:
        if P<=140: rows=[['3/8" (S16)','≤140 bar'],['3/8" (S9)','>140 bar']];rec='3/8" (S16)'
        else:      rows=[['3/8" (S9)','>140 bar']];rec='3/8" (S9)'
    elif D<=21.0:
        if P<=20: rows=[['3/4" (S15)','≤20 bar'],['1" (S14)','>20 bar']];rec='3/4" (S15)'
        else:     rows=[['1" (S14)','>20 bar']];rec='1" (S14)'
    else: rows=[['Special piping','Outside standard range']];rec='Special piping'
    return rows,rec

# ─────────────────────────────────────────────────────────────────────────────
#  Session state
# ─────────────────────────────────────────────────────────────────────────────
for k,v in dict(result=None,error=False,specs=[],
                gas_sel="N2",calc_sel="Pipe Diameter (mm)",field_vals={}).items():
    if k not in st.session_state: st.session_state[k]=v
if not st.session_state.field_vals:
    st.session_state.field_vals={k:str(v) for k,v in DEFAULTS.items()}

# ─────────────────────────────────────────────────────────────────────────────
#  Handle Calculate callback via query params
# ─────────────────────────────────────────────────────────────────────────────
def fkey(f): return f.replace(" ","_").replace("(","").replace(")","").replace("/","").replace("°","deg")

qp=st.query_params
if qp.get("action")=="calc":
    try:
        gas=qp.get("gas","N2"); calc=qp.get("calc","Pipe Diameter (mm)")
        def gv(f): return float(qp.get(fkey(f),DEFAULTS.get(f,0.0)))
        T,Pi,Po,L,D,Q=gv("Temperature (°C)"),gv("Inlet Pressure (bar)"),gv("Outlet Pressure (bar)"),gv("Pipe Length (m)"),gv("Pipe Diameter (mm)"),gv("Flow Rate (LPM)")
        if   calc=="Pipe Diameter (mm)":    r=calc_diameter(Pi,Po,T,L,Q,gas);sp,rec=pipe_opts(r,max(Pi,Po),gas);msg=f"Required Diameter: {r:.2f} mm  ·  Recommended: {rec}"
        elif calc=="Flow Rate (LPM)":       r=calc_flow(Pi,Po,T,L,D,gas);sp,rec=pipe_opts(D,max(Pi,Po),gas);msg=f"Maximum Flow Rate: {r:.1f} L/min  ·  Recommended: {rec}"
        elif calc=="Pipe Length (m)":       r=calc_length(Pi,Po,T,D,Q,gas);sp,rec=pipe_opts(D,max(Pi,Po),gas);msg=f"Maximum Pipe Length: {r:.1f} m  ·  Recommended: {rec}"
        elif calc=="Inlet Pressure (bar)":  r=calc_inlet(Po,T,L,D,Q,gas);sp,rec=pipe_opts(D,max(r,Po),gas);msg=f"Required Inlet Pressure: {r:.2f} bar  ·  Recommended: {rec}"
        elif calc=="Outlet Pressure (bar)": r=calc_outlet(Pi,T,L,D,Q,gas);sp,rec=pipe_opts(D,Pi,gas);msg=f"Estimated Outlet Pressure: {r:.2f} bar  ·  Recommended: {rec}"
        else: msg="Unsupported type.";sp=[]
        st.session_state.update(result=msg,error=False,specs=sp,gas_sel=gas,calc_sel=calc)
        for f in CALC_FIELDS.get(calc,[]):
            if fkey(f) in qp: st.session_state.field_vals[f]=qp[fkey(f)]
    except Exception as e:
        st.session_state.update(result=f"Error: {e}",error=True,specs=[])
    st.query_params.clear()

# ─────────────────────────────────────────────────────────────────────────────
#  Build HTML strings for the panel
# ─────────────────────────────────────────────────────────────────────────────
gas_opts  = "\n".join(f'<option value="{k}"{"selected" if k==st.session_state.gas_sel else ""}>{GAS_LABEL[k]}</option>' for k in GAS_DATA)
calc_opts = "\n".join(f'<option value="{k}"{"selected" if k==st.session_state.calc_sel else ""}>{k}</option>' for k in CALC_FIELDS)

fields_html="".join(f'<div class="frow"><label>{f}</label><input type="number" id="{fkey(f)}" value="{st.session_state.field_vals.get(f,DEFAULTS.get(f,0))}" step="any"></div>'
                    for f in CALC_FIELDS[st.session_state.calc_sel])

def make_result():
    if not st.session_state.result: return ""
    if st.session_state.error:
        return f'<div class="rbox err">⚠ {st.session_state.result}</div>'
    spec_rows="".join(f"<tr><td>{r[0]}</td><td>{r[1]}</td></tr>" for r in st.session_state.specs)
    spec_tbl=f'<table class="stbl"><thead><tr><th>Pipe Spec</th><th>Range</th></tr></thead><tbody>{spec_rows}</tbody></table>' if st.session_state.specs else ""
    return f'<div class="rbox ok">✔ {st.session_state.result}</div>{spec_tbl}'

result_html=make_result()

# ─────────────────────────────────────────────────────────────────────────────
#  Inject:
#   1. HUD as a true position:fixed full-screen iframe (lives in the real DOM)
#   2. Calculator panel as a position:fixed centered overlay (also in the real DOM)
#   3. CSS that strips all Streamlit chrome and makes .stApp transparent
# ─────────────────────────────────────────────────────────────────────────────
inject_html = f"""
<!-- ═══════════════ FULL-SCREEN HUD IFRAME ═══════════════ -->
{"" if not hud_src else f'''
<iframe id="hud-frame" src="{hud_src}"
  style="position:fixed;inset:0;width:100vw;height:100vh;
         border:none;pointer-events:none;z-index:1;"
  sandbox="allow-scripts">
</iframe>
'''}

<!-- ═══════════════ CALCULATOR PANEL ═══════════════ -->
<div id="calc-overlay">
  <div id="calc-panel">
    <h1>⟁ Gas Flow Calculator</h1>

    <div class="srow">
      <label>Gas type</label>
      <select id="gasSelect" onchange="rebuildFields()">{gas_opts}</select>
    </div>
    <div class="srow">
      <label>Calculation type</label>
      <select id="calcSelect" onchange="rebuildFields()">{calc_opts}</select>
    </div>

    <hr class="hdiv">

    <div id="fields">{fields_html}</div>

    <div class="caption">Friction factor f = 0.02 · fixed</div>
    <button id="calc-btn" onclick="doCalc()">▶  Calculate</button>

    <div id="result-area">{result_html}</div>
  </div>
</div>

<!-- ═══════════════ STYLES ═══════════════ -->
<style>
  /* Strip Streamlit chrome */
  #MainMenu, footer, header,
  [data-testid="stToolbar"],
  [data-testid="stDecoration"],
  [data-testid="stStatusWidget"],
  [data-testid="collapsedControl"] {{ display:none !important; }}

  /* Make Streamlit's own containers fully transparent */
  .stApp,
  .stApp > div,
  section.main,
  section.main > div,
  section.main .block-container,
  [data-testid="stAppViewContainer"],
  [data-testid="stAppViewBlockContainer"] {{
    background: transparent !important;
    padding: 0 !important;
    margin: 0 !important;
    max-width: 100vw !important;
    min-height: unset !important;
  }}

  /* Full-screen HUD overlay */
  #calc-overlay {{
    position: fixed;
    inset: 0;
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 1000;
    pointer-events: none;   /* let clicks pass except on the panel */
  }}

  /* The scrollable panel */
  #calc-panel {{
    pointer-events: all;
    font-family: 'Courier New', monospace;
    width: 440px;
    max-width: 90vw;
    max-height: 80vh;          /* ← never taller than 80% of screen */
    overflow-y: auto;          /* ← internal scroll only */
    overflow-x: hidden;
    background: rgba(0,7,6,0.95);
    border: 1.5px solid #00e5cc;
    border-radius: 5px;
    padding: 20px 24px 18px;
    box-shadow: 0 0 30px rgba(0,229,204,0.55), inset 0 0 50px rgba(0,229,204,0.03);
    animation: pglow 2.5s ease-in-out infinite alternate;
    scrollbar-width: thin;
    scrollbar-color: #00e5cc #001a18;
  }}
  #calc-panel::-webkit-scrollbar {{ width: 4px; }}
  #calc-panel::-webkit-scrollbar-track {{ background: #001a18; }}
  #calc-panel::-webkit-scrollbar-thumb {{ background: #00e5cc; border-radius: 2px; }}

  @keyframes pglow {{
    from {{ box-shadow: 0 0 10px rgba(0,229,204,0.3); }}
    to   {{ box-shadow: 0 0 34px rgba(0,229,204,0.7); }}
  }}

  #calc-panel h1 {{
    color: #00ffee;
    text-shadow: 0 0 12px #00e5cc;
    font-size: 12.5px;
    letter-spacing: 3px;
    text-align: center;
    text-transform: uppercase;
    margin-bottom: 14px;
  }}

  .srow {{ display:flex; flex-direction:column; gap:3px; margin-bottom:9px; }}
  .srow label {{ font-size:9.5px; letter-spacing:1.5px; text-transform:uppercase; color:rgba(0,229,204,0.65); }}
  .srow select {{
    background: rgba(0,14,12,0.98); color: #00ffee;
    border: 1px solid rgba(0,229,204,0.38); border-radius: 3px;
    padding: 5px 7px; font-family:'Courier New',monospace; font-size:11px;
    width: 100%; cursor: pointer; outline: none;
  }}
  .srow select:focus {{ border-color:#00e5cc; box-shadow:0 0 7px rgba(0,229,204,0.35); }}

  .hdiv {{ border:none; border-top:1px solid rgba(0,229,204,0.15); margin:10px 0; }}

  .frow {{ display:flex; align-items:center; gap:8px; margin-bottom:7px; }}
  .frow label {{ color:#00e5cc; font-size:11px; width:170px; flex-shrink:0; }}
  .frow input {{
    flex:1; min-width:0;
    background:rgba(0,14,12,0.98); color:#00ffee;
    border:1px solid rgba(0,229,204,0.35); border-radius:3px;
    padding:5px 7px; font-family:'Courier New',monospace; font-size:11px; outline:none;
  }}
  .frow input:focus {{ border-color:#00e5cc; box-shadow:0 0 6px rgba(0,229,204,0.35); }}

  .caption {{ color:rgba(0,229,204,0.4); font-size:9px; text-align:center; margin:3px 0 11px; letter-spacing:1px; }}

  #calc-btn {{
    display:block; width:100%;
    background:rgba(0,229,204,0.07); color:#00ffee;
    border:1.5px solid #00e5cc; border-radius:3px;
    padding:8px 0; font-family:'Courier New',monospace;
    font-size:12px; letter-spacing:3px; cursor:pointer;
    text-transform:uppercase; transition:background .2s, box-shadow .2s;
  }}
  #calc-btn:hover  {{ background:rgba(0,229,204,0.18); box-shadow:0 0 16px rgba(0,229,204,0.5); }}
  #calc-btn:active {{ background:rgba(0,229,204,0.28); }}

  .rbox {{ margin-top:12px; padding:9px 12px; border-radius:3px; font-size:11px; line-height:1.8; }}
  .rbox.ok  {{ background:rgba(0,26,22,0.98); border-left:3px solid #00e5cc; color:#00ffee; }}
  .rbox.err {{ background:rgba(26,0,0,0.98);  border-left:3px solid #ff4444; color:#ff9090; }}

  .stbl {{ width:100%; border-collapse:collapse; margin-top:9px; font-size:10.5px; }}
  .stbl th {{ color:#00e5cc; background:rgba(0,229,204,0.07); border-bottom:1px solid rgba(0,229,204,0.25); padding:4px 8px; text-align:left; letter-spacing:1px; font-size:10px; }}
  .stbl td {{ color:rgba(0,229,204,0.8); border-bottom:1px solid rgba(0,229,204,0.1); padding:4px 8px; }}
</style>

<!-- ═══════════════ CALCULATOR SCRIPT ═══════════════ -->
<script>
const CF = {{
  "Pipe Diameter (mm)":    ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  "Flow Rate (LPM)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  "Pipe Length (m)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Inlet Pressure (bar)":  ["Temperature (°C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Outlet Pressure (bar)": ["Temperature (°C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
}};
const DEF = {{
  "Temperature (°C)":25,"Inlet Pressure (bar)":100,
  "Outlet Pressure (bar)":10,"Pipe Length (m)":10,
  "Pipe Diameter (mm)":10,"Flow Rate (LPM)":100
}};

function toKey(f){{ return f.replace(/ /g,'_').replace(/[()]/g,'').replace(/\\//g,'').replace(/°/g,'deg'); }}

function rebuildFields(){{
  const calc = document.getElementById('calcSelect').value;
  document.getElementById('fields').innerHTML = (CF[calc]||[]).map(f => {{
    const k=toKey(f);
    return `<div class="frow"><label>${{f}}</label><input type="number" id="${{k}}" value="${{DEF[f]??0}}" step="any"></div>`;
  }}).join('');
}}

function doCalc(){{
  const gas=document.getElementById('gasSelect').value;
  const calc=document.getElementById('calcSelect').value;
  const p=new URLSearchParams();
  p.set('action','calc'); p.set('gas',gas); p.set('calc',calc);
  (CF[calc]||[]).forEach(f=>{{
    const el=document.getElementById(toKey(f));
    p.set(toKey(f), el?el.value:(DEF[f]??0));
  }});
  window.parent.location.search='?'+p.toString();
}}
</script>
"""

st.markdown(inject_html, unsafe_allow_html=True)
