"""
High-Pressure Gas Flow Calculator
──────────────────────────────────
• hud_background.html  → full-screen animated background (position:fixed, z-index:1)
• Calculator panel     → fixed centered black box, scrollable, z-index:9999
• The page itself never scrolls – only the panel scrolls internally
• Works with: streamlit run app.py
  (place hud_background.html in the same folder)
"""

import math
import base64
import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Gas Flow Calculator",
    page_icon="💨",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ─────────────────────────────────────────────────────────────────────────────
# DOMAIN DATA
# ─────────────────────────────────────────────────────────────────────────────
GAS_DATA = {
    "N2": 0.028013, "O2": 0.031999, "Ar": 0.039948, "CO2": 0.04401,
    "He": 0.0040026, "H2": 0.002016, "CH4": 0.01604, "C2H2": 0.02604,
    "Forming Gas 1": 0.03881, "Forming Gas 2": 0.02671, "Air": 0.02897,
}
GAS_LABEL = {
    "N2": "N2 – Nitrogen", "O2": "O2 – Oxygen", "Ar": "Ar – Argon",
    "CO2": "CO2 – Carbon Dioxide", "He": "He – Helium", "H2": "H2 – Hydrogen",
    "CH4": "CH4 – Methane", "C2H2": "C2H2 – Acetylene",
    "Forming Gas 1": "H2 3% + Ar 97%  (FG1)",
    "Forming Gas 2": "H2 5% + N2 95%  (FG2)",
    "Air": "Air – Dry Air",
}
CALC_FIELDS = {
    "Pipe Diameter (mm)":    ["Temperature (°C)", "Inlet Pressure (bar)", "Outlet Pressure (bar)", "Pipe Length (m)", "Flow Rate (LPM)"],
    "Flow Rate (LPM)":       ["Temperature (°C)", "Inlet Pressure (bar)", "Outlet Pressure (bar)", "Pipe Length (m)", "Pipe Diameter (mm)"],
    "Pipe Length (m)":       ["Temperature (°C)", "Inlet Pressure (bar)", "Outlet Pressure (bar)", "Pipe Diameter (mm)", "Flow Rate (LPM)"],
    "Inlet Pressure (bar)":  ["Temperature (°C)", "Outlet Pressure (bar)", "Pipe Length (m)", "Pipe Diameter (mm)", "Flow Rate (LPM)"],
    "Outlet Pressure (bar)": ["Temperature (°C)", "Inlet Pressure (bar)", "Pipe Length (m)", "Pipe Diameter (mm)", "Flow Rate (LPM)"],
}
DEFAULTS = {
    "Temperature (°C)": 25.0, "Inlet Pressure (bar)": 100.0,
    "Outlet Pressure (bar)": 10.0, "Pipe Length (m)": 10.0,
    "Pipe Diameter (mm)": 10.0, "Flow Rate (LPM)": 100.0,
}

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICS
# ─────────────────────────────────────────────────────────────────────────────
def _rho(P, T, g):        return (P * GAS_DATA[g]) / (8.314 * T)
def _ra(Pi, Po, T, g):    return (_rho(Pi*1e5, T, g) + _rho(Po*1e5, T, g)) / 2
F = 0.02   # friction factor

def calc_diameter(Pi, Po, Tc, L, Q, g):
    T = Tc + 273.15; Qs = Q / 6e4; dP = (Pi - Po) * 1e5
    if dP <= 0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return ((F * L * 8 * _ra(Pi, Po, T, g) * Qs**2) / (math.pi**2 * dP))**0.2 * 1000

def calc_flow(Pi, Po, Tc, L, D, g):
    T = Tc + 273.15; Dm = D / 1000; dP = (Pi - Po) * 1e5
    if dP <= 0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return math.sqrt((dP * math.pi**2 * Dm**5) / (8 * F * L * _ra(Pi, Po, T, g))) * 6e4

def calc_length(Pi, Po, Tc, D, Q, g):
    T = Tc + 273.15; Dm = D / 1000; Qs = Q / 6e4; dP = (Pi - Po) * 1e5
    if dP <= 0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return (dP * math.pi**2 * Dm**5) / (8 * F * _ra(Pi, Po, T, g) * Qs**2)

def calc_outlet(Pi, Tc, L, D, Q, g):
    T = Tc + 273.15; Dm = D / 1000; Qs = Q / 6e4; Pg = Pi
    for _ in range(20):
        ra = (_rho(Pi*1e5, T, g) + _rho(Pg*1e5, T, g)) / 2
        dP = (8 * F * L * ra * Qs**2) / (math.pi**2 * Dm**5)
        Pn = Pi - dP / 1e5
        if abs(Pn - Pg) < 0.001: return max(Pn, 0.0)
        Pg = Pn
    return max(Pg, 0.0)

def calc_inlet(Po, Tc, L, D, Q, g):
    lo, hi = Po, Po + 100
    while hi < Po + 2000:
        if calc_outlet(hi, Tc, L, D, Q, g) <= Po: break
        hi += 100
    for _ in range(50):
        mid = (lo + hi) / 2; vm = calc_outlet(mid, Tc, L, D, Q, g)
        if abs(vm - Po) < 0.01: return mid
        lo, hi = (mid, hi) if vm > Po else (lo, mid)
    return (lo + hi) / 2

def pipe_opts(D, P, gas):
    if gas == "O2": return [['1" (S14)', 'Required for O2']], '1" (S14)'
    rows, rec = [], None
    if D <= 4.0:
        if P <= 200:    rows = [['1/4" (S6)', '≤200 bar'], ['1/4" (S9)', '≤1379 bar'], ['1/4" (S12)', '>1379 bar']]; rec = '1/4" (S6)'
        elif P <= 1379: rows = [['1/4" (S9)', '≤1379 bar'], ['1/4" (S12)', '>1379 bar']]; rec = '1/4" (S9)'
        else:           rows = [['1/4" (S12)', '>1379 bar']]; rec = '1/4" (S12)'
    elif D <= 7.0:
        if P <= 140: rows = [['3/8" (S16)', '≤140 bar'], ['3/8" (S9)', '>140 bar']]; rec = '3/8" (S16)'
        else:        rows = [['3/8" (S9)', '>140 bar']]; rec = '3/8" (S9)'
    elif D <= 21.0:
        if P <= 20: rows = [['3/4" (S15)', '≤20 bar'], ['1" (S14)', '>20 bar']]; rec = '3/4" (S15)'
        else:       rows = [['1" (S14)', '>20 bar']]; rec = '1" (S14)'
    else: rows = [['Special piping', 'Outside standard range']]; rec = 'Special piping'
    return rows, rec

# ─────────────────────────────────────────────────────────────────────────────
# SESSION STATE
# ─────────────────────────────────────────────────────────────────────────────
for k, v in dict(result=None, error=False, specs=[],
                 gas_sel="N2", calc_sel="Pipe Diameter (mm)", field_vals={}).items():
    if k not in st.session_state:
        st.session_state[k] = v
if not st.session_state.field_vals:
    st.session_state.field_vals = {k: str(v) for k, v in DEFAULTS.items()}

# ─────────────────────────────────────────────────────────────────────────────
# QUERY-PARAM CALLBACK  (Calculate button → parent URL → Streamlit rerun)
# ─────────────────────────────────────────────────────────────────────────────
def fkey(f):
    return f.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "").replace("°", "deg")

qp = st.query_params
if qp.get("action") == "calc":
    try:
        gas  = qp.get("gas",  "N2")
        calc = qp.get("calc", "Pipe Diameter (mm)")
        def gv(f): return float(qp.get(fkey(f), DEFAULTS.get(f, 0.0)))
        T  = gv("Temperature (°C)")
        Pi = gv("Inlet Pressure (bar)")
        Po = gv("Outlet Pressure (bar)")
        L  = gv("Pipe Length (m)")
        D  = gv("Pipe Diameter (mm)")
        Q  = gv("Flow Rate (LPM)")

        if   calc == "Pipe Diameter (mm)":    r = calc_diameter(Pi,Po,T,L,Q,gas); sp,rec = pipe_opts(r,max(Pi,Po),gas); msg = f"Required Diameter: {r:.2f} mm  ·  Recommended: {rec}"
        elif calc == "Flow Rate (LPM)":       r = calc_flow(Pi,Po,T,L,D,gas);     sp,rec = pipe_opts(D,max(Pi,Po),gas); msg = f"Maximum Flow Rate: {r:.1f} L/min  ·  Recommended: {rec}"
        elif calc == "Pipe Length (m)":       r = calc_length(Pi,Po,T,D,Q,gas);   sp,rec = pipe_opts(D,max(Pi,Po),gas); msg = f"Maximum Pipe Length: {r:.1f} m  ·  Recommended: {rec}"
        elif calc == "Inlet Pressure (bar)":  r = calc_inlet(Po,T,L,D,Q,gas);     sp,rec = pipe_opts(D,max(r,Po),gas);  msg = f"Required Inlet Pressure: {r:.2f} bar  ·  Recommended: {rec}"
        elif calc == "Outlet Pressure (bar)": r = calc_outlet(Pi,T,L,D,Q,gas);    sp,rec = pipe_opts(D,Pi,gas);         msg = f"Estimated Outlet Pressure: {r:.2f} bar  ·  Recommended: {rec}"
        else: msg = "Unsupported type."; sp = []

        st.session_state.update(result=msg, error=False, specs=sp, gas_sel=gas, calc_sel=calc)
        for f in CALC_FIELDS.get(calc, []):
            if fkey(f) in qp: st.session_state.field_vals[f] = qp[fkey(f)]
    except Exception as e:
        st.session_state.update(result=f"Error: {e}", error=True, specs=[])
    st.query_params.clear()

# ─────────────────────────────────────────────────────────────────────────────
# LOAD HUD → base64 data-URI
# ─────────────────────────────────────────────────────────────────────────────
def load_hud() -> str:
    p = Path(__file__).resolve().parent / "hud_background.html"
    if not p.exists():
        return ""
    return "data:text/html;base64," + base64.b64encode(p.read_bytes()).decode()

hud_src = load_hud()

# ─────────────────────────────────────────────────────────────────────────────
# BUILD PANEL HTML FRAGMENTS
# ─────────────────────────────────────────────────────────────────────────────
gas_opts  = "\n".join(
    f'<option value="{k}"{"selected" if k==st.session_state.gas_sel else ""}>{GAS_LABEL[k]}</option>'
    for k in GAS_DATA)

calc_opts = "\n".join(
    f'<option value="{k}"{"selected" if k==st.session_state.calc_sel else ""}>{k}</option>'
    for k in CALC_FIELDS)

fields_html = "".join(
    f'<div class="frow"><label>{f}</label>'
    f'<input type="number" id="{fkey(f)}" '
    f'value="{st.session_state.field_vals.get(f, DEFAULTS.get(f, 0))}" step="any"></div>'
    for f in CALC_FIELDS[st.session_state.calc_sel])

def make_result_html():
    if not st.session_state.result:
        return ""
    if st.session_state.error:
        return f'<div class="rbox err">⚠&nbsp;{st.session_state.result}</div>'
    spec_rows = "".join(f"<tr><td>{r[0]}</td><td>{r[1]}</td></tr>" for r in st.session_state.specs)
    spec_tbl  = (f'<table class="stbl"><thead><tr><th>Pipe Spec</th><th>Pressure Range</th></tr>'
                 f'</thead><tbody>{spec_rows}</tbody></table>') if st.session_state.specs else ""
    return f'<div class="rbox ok">✔&nbsp;{st.session_state.result}</div>{spec_tbl}'

result_html = make_result_html()

# ─────────────────────────────────────────────────────────────────────────────
# FULL SELF-CONTAINED HTML  (rendered via components.html at 100 % viewport)
# ─────────────────────────────────────────────────────────────────────────────
# We use components.html with height=0 and a <style> trick to make the iframe
# expand to cover the full Streamlit viewport, then position everything inside.

page_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
/* ── reset ─────────────────────────────────────────── */
*, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
html, body {{
  width: 100%; height: 100%;
  overflow: hidden;          /* page must NEVER scroll */
  background: #000d0d;
  font-family: 'Courier New', monospace;
}}

/* ── HUD layer ─────────────────────────────────────── */
#hud {{
  position: fixed;
  inset: 0;
  width: 100%;
  height: 100%;
  border: none;
  pointer-events: none;
  z-index: 1;
  display: block;
}}

/* ── centered overlay wrapper ──────────────────────── */
#overlay {{
  position: fixed;
  inset: 0;
  display: flex;
  align-items: center;
  justify-content: center;
  z-index: 9999;             /* always above HUD */
  pointer-events: none;      /* pass-through except on #panel */
}}

/* ── the calculator panel ───────────────────────────── */
#panel {{
  pointer-events: all;       /* catch all mouse/touch events */

  /* sizing: fits on screen, never overflows */
  width: 460px;
  max-width: 92vw;
  max-height: 84vh;          /* hard cap – scroll inside if needed */

  /* scrolling: ONLY inside the panel */
  overflow-y: auto;
  overflow-x: hidden;

  /* appearance */
  background: #000000;
  border: 1.5px solid #00e5cc;
  border-radius: 6px;
  padding: 22px 26px 20px;
  box-shadow:
    0 0 0 1px rgba(0,229,204,0.15),
    0 0 30px rgba(0,229,204,0.5),
    inset 0 0 60px rgba(0,229,204,0.03);
  animation: pglow 2.5s ease-in-out infinite alternate;

  /* styled scrollbar */
  scrollbar-width: thin;
  scrollbar-color: #00e5cc #001a18;
}}
#panel::-webkit-scrollbar       {{ width: 5px; }}
#panel::-webkit-scrollbar-track {{ background: #001a18; }}
#panel::-webkit-scrollbar-thumb {{ background: #00e5cc; border-radius: 3px; }}

@keyframes pglow {{
  from {{ box-shadow: 0 0 12px rgba(0,229,204,0.35), inset 0 0 40px rgba(0,229,204,0.02); }}
  to   {{ box-shadow: 0 0 36px rgba(0,229,204,0.7),  inset 0 0 60px rgba(0,229,204,0.05); }}
}}

/* ── panel contents ─────────────────────────────────── */
h1 {{
  color: #00ffee;
  text-shadow: 0 0 14px #00e5cc;
  font-size: 13px;
  letter-spacing: 3px;
  text-align: center;
  text-transform: uppercase;
  margin-bottom: 16px;
  line-height: 1.5;
}}

/* selects */
.srow {{ display: flex; flex-direction: column; gap: 4px; margin-bottom: 10px; }}
.srow label {{
  font-size: 9px; letter-spacing: 2px; text-transform: uppercase;
  color: rgba(0,229,204,0.6);
}}
.srow select {{
  background: #000000; color: #00ffee;
  border: 1px solid rgba(0,229,204,0.4); border-radius: 3px;
  padding: 6px 8px; font-family: 'Courier New', monospace; font-size: 11.5px;
  width: 100%; cursor: pointer; outline: none; appearance: none;
  background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%2300e5cc'/%3E%3C/svg%3E");
  background-repeat: no-repeat; background-position: right 10px center;
  padding-right: 28px;
}}
.srow select:focus {{ border-color: #00e5cc; box-shadow: 0 0 8px rgba(0,229,204,0.4); }}
.srow select option {{ background: #000d0d; color: #00e5cc; }}

/* divider */
.hdiv {{ border: none; border-top: 1px solid rgba(0,229,204,0.18); margin: 12px 0; }}

/* input fields */
.frow {{ display: flex; align-items: center; gap: 10px; margin-bottom: 9px; }}
.frow label {{
  color: #00e5cc; font-size: 11px; letter-spacing: 0.3px;
  width: 180px; flex-shrink: 0; line-height: 1.3;
}}
.frow input {{
  flex: 1; min-width: 0;
  background: #000000; color: #00ffee;
  border: 1px solid rgba(0,229,204,0.35); border-radius: 3px;
  padding: 6px 8px; font-family: 'Courier New', monospace; font-size: 11.5px;
  outline: none;
}}
.frow input:focus {{ border-color: #00e5cc; box-shadow: 0 0 7px rgba(0,229,204,0.4); }}
/* hide number spin buttons */
.frow input::-webkit-inner-spin-button,
.frow input::-webkit-outer-spin-button {{ -webkit-appearance: none; margin: 0; }}
.frow input[type=number] {{ -moz-appearance: textfield; }}

/* caption */
.caption {{
  color: rgba(0,229,204,0.4); font-size: 9px;
  text-align: center; margin: 4px 0 13px; letter-spacing: 1.5px;
}}

/* button */
#calc-btn {{
  display: block; width: 100%;
  background: rgba(0,229,204,0.06); color: #00ffee;
  border: 1.5px solid #00e5cc; border-radius: 4px;
  padding: 9px 0; font-family: 'Courier New', monospace;
  font-size: 12px; letter-spacing: 4px; cursor: pointer;
  text-transform: uppercase;
  transition: background 0.2s, box-shadow 0.2s;
}}
#calc-btn:hover  {{ background: rgba(0,229,204,0.18); box-shadow: 0 0 18px rgba(0,229,204,0.55); }}
#calc-btn:active {{ background: rgba(0,229,204,0.3); }}

/* result */
.rbox {{
  margin-top: 14px; padding: 10px 13px;
  border-radius: 4px; font-size: 11.5px; line-height: 1.9;
  word-break: break-word;
}}
.rbox.ok  {{ background: #000000; border: 1px solid rgba(0,229,204,0.3); border-left: 3px solid #00e5cc; color: #00ffee; }}
.rbox.err {{ background: #000000; border: 1px solid rgba(255,60,60,0.3); border-left: 3px solid #ff4444; color: #ff9090; }}

/* spec table */
.stbl {{ width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 10.5px; }}
.stbl th {{
  color: #00e5cc; background: rgba(0,229,204,0.08);
  border-bottom: 1px solid rgba(0,229,204,0.3);
  padding: 5px 9px; text-align: left; letter-spacing: 1px; font-size: 10px;
  text-transform: uppercase;
}}
.stbl td {{
  color: rgba(0,229,204,0.85);
  border-bottom: 1px solid rgba(0,229,204,0.1);
  padding: 5px 9px;
}}
.stbl tr:last-child td {{ border-bottom: none; }}
</style>
</head>
<body>

<!-- ── HUD background (full screen, behind everything) ── -->
{"<iframe id='hud' src='" + hud_src + "' sandbox='allow-scripts'></iframe>" if hud_src else "<div style='position:fixed;inset:0;background:#000d0d;z-index:1'></div>"}

<!-- ── Centered scrollable calculator panel ── -->
<div id="overlay">
  <div id="panel">

    <h1>⟁ High-Pressure<br>Gas Flow Calculator</h1>

    <div class="srow">
      <label>Gas Type</label>
      <select id="gasSelect" onchange="rebuildFields()">
        {gas_opts}
      </select>
    </div>

    <div class="srow">
      <label>Calculation Type</label>
      <select id="calcSelect" onchange="rebuildFields()">
        {calc_opts}
      </select>
    </div>

    <hr class="hdiv">

    <div id="fields">
      {fields_html}
    </div>

    <div class="caption">Friction factor &nbsp;f = 0.02 &nbsp;·&nbsp; fixed constant</div>

    <button id="calc-btn" onclick="doCalc()">▶ &nbsp; Calculate</button>

    <div id="result-area">{result_html}</div>

  </div><!-- #panel -->
</div><!-- #overlay -->

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

function toKey(f){{
  return f.replace(/ /g,'_').replace(/[()]/g,'').replace(/\\//g,'').replace(/°/g,'deg');
}}

function rebuildFields(){{
  const calc = document.getElementById('calcSelect').value;
  document.getElementById('fields').innerHTML = (CF[calc]||[]).map(f => {{
    const k = toKey(f);
    return `<div class="frow">
      <label>${{f}}</label>
      <input type="number" id="${{k}}" value="${{DEF[f]??0}}" step="any">
    </div>`;
  }}).join('');
}}

function doCalc(){{
  const gas  = document.getElementById('gasSelect').value;
  const calc = document.getElementById('calcSelect').value;
  const p    = new URLSearchParams();
  p.set('action','calc');
  p.set('gas', gas);
  p.set('calc', calc);
  (CF[calc]||[]).forEach(f => {{
    const el = document.getElementById(toKey(f));
    p.set(toKey(f), el ? el.value : (DEF[f]??0));
  }});
  // Trigger Streamlit rerun via parent window navigation
  window.parent.location.search = '?' + p.toString();
}}
</script>
</body>
</html>"""

# ─────────────────────────────────────────────────────────────────────────────
# RENDER
# We set the iframe to cover 100 % of the Streamlit viewport by injecting
# CSS into the parent page as well.
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
  /* strip ALL Streamlit chrome */
  #MainMenu, footer, header,
  [data-testid="stToolbar"],
  [data-testid="stDecoration"],
  [data-testid="stStatusWidget"],
  [data-testid="collapsedControl"] { display: none !important; }

  /* make the Streamlit app itself invisible (our iframe fills everything) */
  .stApp,
  [data-testid="stAppViewContainer"],
  section.main {
    background: #000d0d !important;
    padding: 0 !important;
    margin: 0 !important;
  }
  section.main .block-container {
    padding: 0 !important;
    margin: 0 !important;
    max-width: 100vw !important;
  }

  /* expand the components.html iframe to fill the whole viewport */
  iframe[title="streamlit_components_v1.html_v1"] {
    position: fixed !important;
    inset: 0 !important;
    width: 100vw !important;
    height: 100vh !important;
    border: none !important;
    z-index: 100 !important;
  }
</style>
""", unsafe_allow_html=True)

# Render page_html inside a components iframe that covers 100 % of the screen
components.html(page_html, height=800, scrolling=False)
