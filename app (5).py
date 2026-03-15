import math
import re
import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="High-Pressure Gas Flow Calculator",
    page_icon="💨",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Hide ALL Streamlit chrome ─────────────────────────────────────────────────
st.markdown("""
<style>
  #MainMenu, footer, header,
  [data-testid="stToolbar"],
  [data-testid="stDecoration"],
  [data-testid="stStatusWidget"] { display:none !important; }
  .stApp { background: #000d0d !important; }
  section.main .block-container {
    padding: 0 !important;
    margin: 0 !important;
    max-width: 100vw !important;
    width: 100vw !important;
  }
</style>
""", unsafe_allow_html=True)

# ── Load HUD HTML ─────────────────────────────────────────────────────────────
def load_hud_html() -> str:
    path = Path(__file__).resolve().parent / "hud_background.html"
    return path.read_text(encoding="utf-8") if path.exists() else ""

hud_raw = load_hud_html()

# ── Data ──────────────────────────────────────────────────────────────────────
GAS_DATA = {
    "N2": 0.028013, "O2": 0.031999, "Ar": 0.039948, "CO2": 0.04401,
    "He": 0.0040026, "H2": 0.002016, "CH4": 0.01604, "C2H2": 0.02604,
    "Forming Gas 1": 0.03881, "Forming Gas 2": 0.02671, "Air": 0.02897,
}
GAS_DISPLAY = {
    "N2": "N2 – Nitrogen", "O2": "O2 – Oxygen", "Ar": "Ar – Argon",
    "CO2": "CO2 – Carbon Dioxide", "He": "He – Helium", "H2": "H2 – Hydrogen",
    "CH4": "CH4 – Methane", "C2H2": "C2H2 – Acetylene",
    "Forming Gas 1": "H2-3% + Ar-97%  (FG1)",
    "Forming Gas 2": "H2-5% + N2-95%  (FG2)",
    "Air": "Air – Dry Air",
}
CALC_FIELDS = {
    "Pipe Diameter (mm)":    ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
    "Flow Rate (LPM)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
    "Pipe Length (m)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
    "Inlet Pressure (bar)":  ["Temperature (°C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
    "Outlet Pressure (bar)": ["Temperature (°C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
}
DEFAULTS = {
    "Temperature (°C)": 25.0, "Inlet Pressure (bar)": 100.0,
    "Outlet Pressure (bar)": 10.0, "Pipe Length (m)": 10.0,
    "Pipe Diameter (mm)": 10.0, "Flow Rate (LPM)": 100.0,
}
FRICTION_FACTOR = 0.02

# ── Physics ───────────────────────────────────────────────────────────────────
def rho(P_pa, T_k, gas): return (P_pa * GAS_DATA[gas]) / (8.314 * T_k)
def rho_avg(Pi, Po, Tk, g): return (rho(Pi*1e5,Tk,g)+rho(Po*1e5,Tk,g))/2

def calc_diameter(Pi,Po,Tc,L,Q,g):
    Tk=Tc+273.15;Qs=Q/6e4;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet must exceed outlet pressure.")
    return ((FRICTION_FACTOR*L*8*rho_avg(Pi,Po,Tk,g)*Qs**2)/(math.pi**2*dP))**0.2*1000

def calc_flow(Pi,Po,Tc,L,D,g):
    Tk=Tc+273.15;Dm=D/1000;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet must exceed outlet pressure.")
    return math.sqrt((dP*math.pi**2*Dm**5)/(8*FRICTION_FACTOR*L*rho_avg(Pi,Po,Tk,g)))*6e4

def calc_length(Pi,Po,Tc,D,Q,g):
    Tk=Tc+273.15;Dm=D/1000;Qs=Q/6e4;dP=(Pi-Po)*1e5
    if dP<=0: raise ValueError("Inlet must exceed outlet pressure.")
    return (dP*math.pi**2*Dm**5)/(8*FRICTION_FACTOR*rho_avg(Pi,Po,Tk,g)*Qs**2)

def calc_outlet(Pi,Tc,L,D,Q,g):
    Tk=Tc+273.15;Dm=D/1000;Qs=Q/6e4;Pg=Pi
    for _ in range(20):
        ra=(rho(Pi*1e5,Tk,g)+rho(Pg*1e5,Tk,g))/2
        dP=(8*FRICTION_FACTOR*L*ra*Qs**2)/(math.pi**2*Dm**5)
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
        mid=(lo+hi)/2
        vm=calc_outlet(mid,Tc,L,D,Q,g)
        if abs(vm-Po)<0.01: return mid
        if vm>Po: lo=mid
        else: hi=mid
    return (lo+hi)/2

def pipe_options(D,P,gas):
    if gas=="O2": return [["1\" (Spec S14)","Required for O2"]],"1\" (Spec S14)"
    rows,rec=[],None
    if D<=4.0:
        if P<=200:   rows=[['1/4" (S6)','≤200 bar'],['1/4" (S9)','≤1379 bar'],['1/4" (S12)','>1379 bar']];rec='1/4" (S6)'
        elif P<=1379:rows=[['1/4" (S9)','≤1379 bar'],['1/4" (S12)','>1379 bar']];rec='1/4" (S9)'
        else:        rows=[['1/4" (S12)','>1379 bar']];rec='1/4" (S12)'
    elif D<=7.0:
        if P<=140:   rows=[['3/8" (S16)','≤140 bar'],['3/8" (S9)','>140 bar']];rec='3/8" (S16)'
        else:        rows=[['3/8" (S9)','>140 bar']];rec='3/8" (S9)'
    elif D<=21.0:
        if P<=20:    rows=[['3/4" (S15)','≤20 bar'],['1" (S14)','>20 bar']];rec='3/4" (S15)'
        else:        rows=[['1" (S14)','>20 bar']];rec='1" (S14)'
    else:            rows=[['Special piping','Outside standard range']];rec='Special piping'
    return rows,rec

# ── Session state ─────────────────────────────────────────────────────────────
for k,v in [("result",None),("error",False),("specs",[]),
            ("gas_sel","N2"),("calc_sel","Pipe Diameter (mm)"),("field_vals",{})]:
    if k not in st.session_state: st.session_state[k] = v
if not st.session_state.field_vals:
    st.session_state.field_vals = {k: str(v) for k,v in DEFAULTS.items()}

# ── Handle query-param callback from the iframe ───────────────────────────────
qp = st.query_params
if qp.get("action") == "calc":
    try:
        gas  = qp.get("gas","N2")
        calc = qp.get("calc","Pipe Diameter (mm)")
        def gv(f): 
            key=f.replace(" ","_").replace("(","").replace(")","").replace("/","").replace("°","deg")
            return float(qp.get(key, DEFAULTS.get(f,0.0)))
        T,Pi,Po,L,D,Q = gv("Temperature (°C)"),gv("Inlet Pressure (bar)"),gv("Outlet Pressure (bar)"),gv("Pipe Length (m)"),gv("Pipe Diameter (mm)"),gv("Flow Rate (LPM)")
        if   calc=="Pipe Diameter (mm)":    r=calc_diameter(Pi,Po,T,L,Q,gas);  sp,rec=pipe_options(r,max(Pi,Po),gas); msg=f"Required Diameter: {r:.2f} mm  ·  Recommended: {rec}"
        elif calc=="Flow Rate (LPM)":       r=calc_flow(Pi,Po,T,L,D,gas);     sp,rec=pipe_options(D,max(Pi,Po),gas); msg=f"Maximum Flow Rate: {r:.1f} L/min  ·  Recommended: {rec}"
        elif calc=="Pipe Length (m)":       r=calc_length(Pi,Po,T,D,Q,gas);   sp,rec=pipe_options(D,max(Pi,Po),gas); msg=f"Maximum Pipe Length: {r:.1f} m  ·  Recommended: {rec}"
        elif calc=="Inlet Pressure (bar)":  r=calc_inlet(Po,T,L,D,Q,gas);     sp,rec=pipe_options(D,max(r,Po),gas);  msg=f"Required Inlet Pressure: {r:.2f} bar  ·  Recommended: {rec}"
        elif calc=="Outlet Pressure (bar)": r=calc_outlet(Pi,T,L,D,Q,gas);    sp,rec=pipe_options(D,Pi,gas);         msg=f"Estimated Outlet Pressure: {r:.2f} bar  ·  Recommended: {rec}"
        else: msg="Unsupported type."; sp=[]
        st.session_state.update(result=msg, error=False, specs=sp, gas_sel=gas, calc_sel=calc)
        for f in CALC_FIELDS.get(calc,[]):
            key=f.replace(" ","_").replace("(","").replace(")","").replace("/","").replace("°","deg")
            if key in qp: st.session_state.field_vals[f]=qp[key]
    except Exception as e:
        st.session_state.update(result=f"Error: {e}", error=True, specs=[])
    st.query_params.clear()

# ── Extract parts from HUD HTML ───────────────────────────────────────────────
hud_styles  = "\n".join(re.findall(r'<style[^>]*>(.*?)</style>', hud_raw, re.DOTALL))
hud_body_m  = re.search(r'<body[^>]*>(.*?)</body>', hud_raw, re.DOTALL)
hud_content = hud_body_m.group(1) if hud_body_m else ""
hud_scripts = "\n".join(re.findall(r'<script[^>]*>(.*?)</script>', hud_raw, re.DOTALL))

# ── Build select options ──────────────────────────────────────────────────────
gas_opts  = "\n".join(f'<option value="{k}"{"selected" if k==st.session_state.gas_sel else ""}>{GAS_DISPLAY[k]}</option>' for k in GAS_DATA)
calc_opts = "\n".join(f'<option value="{k}"{"selected" if k==st.session_state.calc_sel else ""}>{k}</option>' for k in CALC_FIELDS)

# ── Input fields HTML ─────────────────────────────────────────────────────────
def field_key(f): return f.replace(" ","_").replace("(","").replace(")","").replace("/","").replace("°","deg")
fields_html = "".join(f"""
  <div class="field-row">
    <label>{f}</label>
    <input type="number" id="{field_key(f)}" value="{st.session_state.field_vals.get(f,DEFAULTS.get(f,0))}" step="any">
  </div>""" for f in CALC_FIELDS[st.session_state.calc_sel])

# ── Result HTML ───────────────────────────────────────────────────────────────
result_html = ""
if st.session_state.result:
    if st.session_state.error:
        result_html = f'<div class="result-box err">⚠ {st.session_state.result}</div>'
    else:
        spec_rows = "".join(f"<tr><td>{r[0]}</td><td>{r[1]}</td></tr>" for r in st.session_state.specs)
        spec_table = f'<table class="spec-table"><thead><tr><th>Pipe Spec</th><th>Range</th></tr></thead><tbody>{spec_rows}</tbody></table>' if st.session_state.specs else ""
        result_html = f'<div class="result-box ok">✔ {st.session_state.result}</div>{spec_table}'

# ── Full HTML page ────────────────────────────────────────────────────────────
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
/* ══ HUD original styles ══ */
{hud_styles}

/* ══ Page reset ══ */
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
html,body{{width:100%;height:100%;background:#000d0d;overflow:hidden;font-family:'Courier New',monospace;color:#00e5cc}}

/* ══ HUD layer (fixed, behind everything) ══ */
.hud-root{{position:fixed;inset:0;width:100vw;height:100vh;pointer-events:none;z-index:1}}

/* ══ Calculator wrapper (centered, above HUD) ══ */
#calc-wrap{{
  position:fixed;inset:0;
  display:flex;align-items:center;justify-content:center;
  z-index:500;
  pointer-events:none;
}}

/* ══ The panel itself ══ */
#calc-panel{{
  pointer-events:all;
  background:rgba(0,6,6,0.94);
  border:1.5px solid #00e5cc;
  border-radius:5px;
  width:min(460px,88vw);
  max-height:88vh;
  overflow-y:auto;
  padding:20px 24px 18px;
  box-shadow:0 0 28px rgba(0,229,204,0.55),inset 0 0 50px rgba(0,229,204,0.03);
  animation:pglow 2.5s ease-in-out infinite alternate;
  scrollbar-width:thin;scrollbar-color:#00e5cc #000;
}}
#calc-panel::-webkit-scrollbar{{width:4px}}
#calc-panel::-webkit-scrollbar-thumb{{background:#00e5cc;border-radius:2px}}
@keyframes pglow{{
  from{{box-shadow:0 0 10px rgba(0,229,204,0.3)}}
  to  {{box-shadow:0 0 32px rgba(0,229,204,0.65)}}
}}

/* ══ Title ══ */
#calc-panel h1{{
  color:#00ffee;text-shadow:0 0 12px #00e5cc;
  font-size:13px;letter-spacing:3px;text-align:center;
  text-transform:uppercase;margin-bottom:14px;
}}

/* ══ Select rows ══ */
.sel-row{{display:flex;flex-direction:column;gap:3px;margin-bottom:9px}}
.sel-row label{{font-size:10px;letter-spacing:1.5px;text-transform:uppercase;color:rgba(0,229,204,0.7)}}
.sel-row select{{
  background:rgba(0,15,13,0.97);color:#00ffee;
  border:1px solid rgba(0,229,204,0.4);border-radius:3px;
  padding:5px 7px;font-family:'Courier New',monospace;font-size:11px;
  width:100%;cursor:pointer;outline:none;
}}
.sel-row select:focus{{border-color:#00e5cc;box-shadow:0 0 7px rgba(0,229,204,0.4)}}

/* ══ Divider ══ */
.div{{border:none;border-top:1px solid rgba(0,229,204,0.18);margin:10px 0}}

/* ══ Input field rows ══ */
.field-row{{display:flex;align-items:center;gap:8px;margin-bottom:7px}}
.field-row label{{
  color:#00e5cc;font-size:11px;letter-spacing:0.3px;
  width:175px;flex-shrink:0;
}}
.field-row input{{
  flex:1;min-width:0;
  background:rgba(0,15,13,0.97);color:#00ffee;
  border:1px solid rgba(0,229,204,0.38);border-radius:3px;
  padding:5px 7px;font-family:'Courier New',monospace;font-size:11px;outline:none;
}}
.field-row input:focus{{border-color:#00e5cc;box-shadow:0 0 6px rgba(0,229,204,0.4)}}

/* ══ Caption ══ */
.caption{{color:rgba(0,229,204,0.45);font-size:9.5px;text-align:center;margin:3px 0 11px;letter-spacing:1px}}

/* ══ Button ══ */
#calc-btn{{
  display:block;width:100%;
  background:rgba(0,229,204,0.07);color:#00ffee;
  border:1.5px solid #00e5cc;border-radius:3px;
  padding:8px 0;font-family:'Courier New',monospace;
  font-size:12px;letter-spacing:3px;cursor:pointer;
  transition:background .2s,box-shadow .2s;text-transform:uppercase;
}}
#calc-btn:hover{{background:rgba(0,229,204,0.18);box-shadow:0 0 16px rgba(0,229,204,0.55)}}
#calc-btn:active{{background:rgba(0,229,204,0.28)}}

/* ══ Result ══ */
.result-box{{margin-top:12px;padding:9px 12px;border-radius:3px;font-size:11.5px;line-height:1.75}}
.result-box.ok {{background:rgba(0,28,24,0.97);border-left:3px solid #00e5cc;color:#00ffee}}
.result-box.err{{background:rgba(28,0,0,0.97);border-left:3px solid #ff4444;color:#ff9090}}

/* ══ Spec table ══ */
.spec-table{{width:100%;border-collapse:collapse;margin-top:9px;font-size:10.5px}}
.spec-table th{{
  color:#00e5cc;background:rgba(0,229,204,0.07);
  border-bottom:1px solid rgba(0,229,204,0.28);
  padding:4px 8px;text-align:left;letter-spacing:1px;font-size:10px;
}}
.spec-table td{{color:rgba(0,229,204,0.82);border-bottom:1px solid rgba(0,229,204,0.1);padding:4px 8px}}
</style>
</head>
<body>

<!-- HUD background -->
<div class="hud-root">
{hud_content}
</div>

<!-- Calculator panel -->
<div id="calc-wrap">
  <div id="calc-panel">
    <h1>⟁ High-Pressure Gas Flow Calculator</h1>

    <div class="sel-row">
      <label>Gas type</label>
      <select id="gasSelect" onchange="updateFields()">
        {gas_opts}
      </select>
    </div>

    <div class="sel-row">
      <label>Calculation type</label>
      <select id="calcSelect" onchange="updateFields()">
        {calc_opts}
      </select>
    </div>

    <hr class="div">

    <div id="inputFields">
      {fields_html}
    </div>

    <div class="caption">Friction factor (f) = 0.02 · fixed constant</div>
    <button id="calc-btn" onclick="submitCalc()">▶ Calculate</button>

    <div id="resultArea">{result_html}</div>
  </div>
</div>

<script>
/* ── HUD scripts ── */
{hud_scripts}

/* ── Calculator logic ── */
const CF = {{
  "Pipe Diameter (mm)":    ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  "Flow Rate (LPM)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  "Pipe Length (m)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Inlet Pressure (bar)":  ["Temperature (°C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Outlet Pressure (bar)": ["Temperature (°C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
}};
const DEF = {{"Temperature (°C)":25,"Inlet Pressure (bar)":100,"Outlet Pressure (bar)":10,"Pipe Length (m)":10,"Pipe Diameter (mm)":10,"Flow Rate (LPM)":100}};

function toKey(f){{return f.replace(/ /g,'_').replace(/[()]/g,'').replace(/\\//g,'').replace(/°/g,'deg')}}

function updateFields(){{
  const calc = document.getElementById('calcSelect').value;
  const fields = CF[calc]||[];
  document.getElementById('inputFields').innerHTML = fields.map(f=>{{
    const k=toKey(f);
    return `<div class="field-row"><label>${{f}}</label><input type="number" id="${{k}}" value="${{DEF[f]??0}}" step="any"></div>`;
  }}).join('');
}}

function submitCalc(){{
  const gas  = document.getElementById('gasSelect').value;
  const calc = document.getElementById('calcSelect').value;
  const fields = CF[calc]||[];
  const p = new URLSearchParams();
  p.set('action','calc'); p.set('gas',gas); p.set('calc',calc);
  fields.forEach(f=>{{
    const el=document.getElementById(toKey(f));
    p.set(toKey(f), el?el.value:DEF[f]??0);
  }});
  window.parent.location.search='?'+p.toString();
}}
</script>
</body>
</html>"""

components.html(html, height=720, scrolling=False)
