"""
High-Pressure Gas Flow Calculator
──────────────────────────────────
Layout:
  • HUD animated background fills entire screen  (z-index 1)
  • LEFT panel  – black, fit-content, scrollable:
      Title "HIGH-PRESSURE GAS FLOW CALCULATOR"  ← top of this panel
      Gas type dropdown
      Calculation type dropdown
      All input fields
  • RIGHT panel – black, fit-content, scrollable:
      Friction-factor caption
      CALCULATE button
      Green bold result box
      Pipe-spec table
  • HUD graphics are fully visible: top strip, bottom strip,
    far-left strip, far-right strip, centre gap between panels
  • Nothing ever overlaps the HUD graphics

Run:  streamlit run app.py
Req:  hud_background.html in the same folder
"""

import math
import base64
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

# Strip every pixel of Streamlit's own UI; expand the components iframe
st.markdown("""
<style>
  #MainMenu, footer, header,
  [data-testid="stToolbar"],
  [data-testid="stDecoration"],
  [data-testid="stStatusWidget"],
  [data-testid="collapsedControl"] { display: none !important; }
  .stApp { background: #000d0d !important; }
  section.main .block-container {
    padding: 0 !important; margin: 0 !important;
    max-width: 100vw !important;
  }
  iframe[title="streamlit_components_v1.html_v1"] {
    position: fixed !important;
    inset: 0 !important;
    width:  100vw !important;
    height: 100vh !important;
    border: none !important;
    z-index: 100 !important;
  }
</style>
""", unsafe_allow_html=True)

# ── Domain data ───────────────────────────────────────────────────────────────
GAS_DATA = {
    "N2": 0.028013, "O2": 0.031999, "Ar": 0.039948, "CO2": 0.04401,
    "He": 0.0040026, "H2": 0.002016, "CH4": 0.01604, "C2H2": 0.02604,
    "Forming Gas 1": 0.03881, "Forming Gas 2": 0.02671, "Air": 0.02897,
}
GAS_DISPLAY = {
    "N2": "N2 (Nitrogen)", "O2": "O2 (Oxygen)", "Ar": "Ar (Argon)",
    "CO2": "CO2 (Carbon Dioxide)", "He": "He (Helium)", "H2": "H2 (Hydrogen)",
    "CH4": "CH4 (Methane)", "C2H2": "C2H2 (Acetylene)",
    "Forming Gas 1": "(H2-3% + Ar-97%) Forming Gas 1",
    "Forming Gas 2": "(H2-5% + N2-95%) Forming Gas 2",
    "Air": "Air (Dry Air)",
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
FRICTION_FACTOR = 0.02

# ── Physics ───────────────────────────────────────────────────────────────────
def igd(P, T, g):
    return (P * GAS_DATA[g]) / (8.314 * T)

def ra(Pi, Po, T, g):
    return (igd(Pi * 1e5, T, g) + igd(Po * 1e5, T, g)) / 2.0

def calc_required_diameter(Pi, Po, Tc, L, Q, g):
    T = Tc + 273.15; Qs = Q / 6e4; dP = (Pi - Po) * 1e5
    if dP <= 0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return ((FRICTION_FACTOR * L * 8.0 * ra(Pi, Po, T, g) * Qs**2) /
            (math.pi**2 * dP)) ** 0.2 * 1000.0

def calc_flow_rate(Pi, Po, Tc, L, D, g):
    T = Tc + 273.15; Dm = D / 1000.0; dP = (Pi - Po) * 1e5
    if dP <= 0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return math.sqrt((dP * math.pi**2 * Dm**5) /
                     (8.0 * FRICTION_FACTOR * L * ra(Pi, Po, T, g))) * 6e4

def calc_max_length(Pi, Po, Tc, D, Q, g):
    T = Tc + 273.15; Dm = D / 1000.0; Qs = Q / 6e4; dP = (Pi - Po) * 1e5
    if dP <= 0: raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return (dP * math.pi**2 * Dm**5) / (8.0 * FRICTION_FACTOR * ra(Pi, Po, T, g) * Qs**2)

def calc_outlet_pressure(Pi, Tc, L, D, Q, g):
    """Bisection — guaranteed convergence, no oscillation."""
    T = Tc + 273.15; Dm = D / 1000.0; Qs = Q / 6e4
    def res(Po):
        return (Pi - Po) * 1e5 - (8.0 * FRICTION_FACTOR * L * ra(Pi, Po, T, g) * Qs**2) / (math.pi**2 * Dm**5)
    lo, hi = 0.0, Pi
    for _ in range(60):
        if abs(hi - lo) < 1e-4: break
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if res(mid) > 0 else (lo, mid)
    return max((lo + hi) / 2, 0.0)

def calc_required_inlet_pressure(Po, Tc, L, D, Q, g):
    """Outlet grows with inlet → expand hi until outlet(hi) ≥ Po, then bisect."""
    lo = Po; hi = Po + 10.0
    while hi < Po + 2000:
        if calc_outlet_pressure(hi, Tc, L, D, Q, g) >= Po: break
        hi += 10.0
    for _ in range(60):
        mid = (lo + hi) / 2
        vm  = calc_outlet_pressure(mid, Tc, L, D, Q, g)
        if abs(vm - Po) < 0.005: return mid
        if vm < Po: lo = mid
        else:       hi = mid
    return (lo + hi) / 2

def determine_pipe_options(D, P, gas):
    if gas == "O2":
        return ['1" tube (Spec S14) – required for O2'], '1" tube (Spec S14)'
    opts, rec = [], None
    if D <= 4.0:
        if P <= 200:
            opts = ['1/4" tube (Spec S6) – up to 200 bar',
                    '1/4" tube (Spec S9) – up to 1379 bar',
                    '1/4" tube (Spec S12) – above 1379 bar']; rec = '1/4" tube (Spec S6)'
        elif P <= 1379:
            opts = ['1/4" tube (Spec S9) – up to 1379 bar',
                    '1/4" tube (Spec S12) – above 1379 bar']; rec = '1/4" tube (Spec S9)'
        else:
            opts = ['1/4" tube (Spec S12) – above 1379 bar']; rec = '1/4" tube (Spec S12)'
    elif D <= 7.0:
        if P <= 140:
            opts = ['3/8" tube (Spec S16) – up to 140 bar',
                    '3/8" tube (Spec S9) – above 140 bar'];   rec = '3/8" tube (Spec S16)'
        else:
            opts = ['3/8" tube (Spec S9) – above 140 bar'];   rec = '3/8" tube (Spec S9)'
    elif D <= 21.0:
        if P <= 20:
            opts = ['3/4" tube (Spec S15) – up to 20 bar',
                    '1" tube (Spec S14) – above 20 bar'];     rec = '3/4" tube (Spec S15)'
        else:
            opts = ['1" tube (Spec S14) – above 20 bar'];     rec = '1" tube (Spec S14)'
    else:
        opts = ['Special piping required (outside standard range)']; rec = 'Special piping'
    return opts, rec

# ── Session state ─────────────────────────────────────────────────────────────
for k, v in dict(result=None, error=False, specs=[],
                 gas_sel="N2", calc_sel="Pipe Diameter (mm)", field_vals={}).items():
    if k not in st.session_state:
        st.session_state[k] = v
if not st.session_state.field_vals:
    st.session_state.field_vals = {k: str(v) for k, v in DEFAULTS.items()}

def fkey(f):
    return (f.replace(" ", "_").replace("(", "").replace(")", "")
             .replace("/", "").replace("°", "deg"))

# ── Query-param callback ──────────────────────────────────────────────────────
qp = st.query_params
if qp.get("action") == "calc":
    try:
        gas  = qp.get("gas",  "N2")
        calc = qp.get("calc", "Pipe Diameter (mm)")
        def gv(f): return float(qp.get(fkey(f), DEFAULTS.get(f, 0.0)))
        Tc = gv("Temperature (°C)");      Pi = gv("Inlet Pressure (bar)")
        Po = gv("Outlet Pressure (bar)"); L  = gv("Pipe Length (m)")
        D  = gv("Pipe Diameter (mm)");    Q  = gv("Flow Rate (LPM)")

        if calc == "Pipe Diameter (mm)":
            r = calc_required_diameter(Pi, Po, Tc, L, Q, gas)
            sp, rec = determine_pipe_options(r, max(Pi, Po), gas)
            msg = f"Required Diameter: {r:.2f} mm<br><small>Recommended: {rec}</small>"
        elif calc == "Flow Rate (LPM)":
            r = calc_flow_rate(Pi, Po, Tc, L, D, gas)
            sp, rec = determine_pipe_options(D, max(Pi, Po), gas)
            msg = f"Maximum Flow Rate: {r:.1f} L/min<br><small>Recommended: {rec}</small>"
        elif calc == "Pipe Length (m)":
            r = calc_max_length(Pi, Po, Tc, D, Q, gas)
            sp, rec = determine_pipe_options(D, max(Pi, Po), gas)
            msg = f"Maximum Pipe Length: {r:.1f} m<br><small>Recommended: {rec}</small>"
        elif calc == "Inlet Pressure (bar)":
            r = calc_required_inlet_pressure(Po, Tc, L, D, Q, gas)
            sp, rec = determine_pipe_options(D, max(r, Po), gas)
            msg = f"Required Inlet Pressure: {r:.2f} bar<br><small>Recommended: {rec}</small>"
        elif calc == "Outlet Pressure (bar)":
            r = calc_outlet_pressure(Pi, Tc, L, D, Q, gas)
            sp, rec = determine_pipe_options(D, Pi, gas)
            msg = f"Estimated Outlet Pressure: {r:.2f} bar<br><small>Recommended: {rec}</small>"
        else:
            msg = "Unsupported calculation type."; sp = []

        st.session_state.update(result=msg, error=False, specs=sp,
                                gas_sel=gas, calc_sel=calc)
        for f in CALC_FIELDS.get(calc, []):
            if fkey(f) in qp:
                st.session_state.field_vals[f] = qp[fkey(f)]
    except Exception as e:
        st.session_state.update(result=f"Error: {e}", error=True, specs=[])
    st.query_params.clear()

# ── Load HUD ──────────────────────────────────────────────────────────────────
def load_hud() -> str:
    p = Path(__file__).resolve().parent / "hud_background.html"
    if not p.exists():
        return ""
    return "data:text/html;base64," + base64.b64encode(p.read_bytes()).decode()

hud_src = load_hud()

# ── HTML fragments ────────────────────────────────────────────────────────────
gas_opts = "\n".join(
    f'<option value="{k}"{"selected" if k == st.session_state.gas_sel else ""}>'
    f'{GAS_DISPLAY[k]}</option>'
    for k in GAS_DATA)

calc_opts = "\n".join(
    f'<option value="{k}"{"selected" if k == st.session_state.calc_sel else ""}>'
    f'{k}</option>'
    for k in CALC_FIELDS)

def render_field(f):
    v = st.session_state.field_vals.get(f, DEFAULTS.get(f, 0))
    k = fkey(f)
    return (f'<div class="frow">'
            f'<label>{f}</label>'
            f'<input type="number" id="{k}" value="{v}" step="any">'
            f'</div>')

fields_html = "".join(render_field(f) for f in CALC_FIELDS[st.session_state.calc_sel])

def make_result_html():
    if not st.session_state.result:
        return ""
    if st.session_state.error:
        return f'<div class="rbox err">⚠&nbsp;{st.session_state.result}</div>'
    rows = ""
    for spec in st.session_state.specs:
        if "–" in spec:
            name, detail = spec.split("–", 1)
            rows += f"<tr><td>{name.strip()}</td><td>{detail.strip()}</td></tr>"
        else:
            rows += f"<tr><td>{spec}</td><td></td></tr>"
    tbl = (f'<p class="sp-title">Possible pipe specifications:</p>'
           f'<table class="stbl">'
           f'<thead><tr><th>Pipe Spec</th><th>Details</th></tr></thead>'
           f'<tbody>{rows}</tbody></table>') if st.session_state.specs else ""
    return f'<div class="rbox ok">{st.session_state.result}</div>{tbl}'

result_html = make_result_html()

# ── JS maps ───────────────────────────────────────────────────────────────────
CF_JS = """{
  "Pipe Diameter (mm)":    ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Flow Rate (LPM)"],
  "Flow Rate (LPM)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)"],
  "Pipe Length (m)":       ["Temperature (°C)","Inlet Pressure (bar)","Outlet Pressure (bar)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Inlet Pressure (bar)":  ["Temperature (°C)","Outlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"],
  "Outlet Pressure (bar)": ["Temperature (°C)","Inlet Pressure (bar)","Pipe Length (m)","Pipe Diameter (mm)","Flow Rate (LPM)"]
}"""
DEF_JS = """{
  "Temperature (°C)":25,"Inlet Pressure (bar)":100,
  "Outlet Pressure (bar)":10,"Pipe Length (m)":10,
  "Pipe Diameter (mm)":10,"Flow Rate (LPM)":100
}"""

# ── Full page HTML ────────────────────────────────────────────────────────────
page_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
/* ══ reset ══ */
*,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
html,body{{
  width:100%;height:100%;
  overflow:hidden;        /* page NEVER scrolls — only panels scroll */
  background:#000d0d;
  font-family:'Courier New',monospace;
  color:#00e5cc;
}}

/* ══ HUD — full-screen, below everything ══ */
#hud{{
  position:fixed;inset:0;
  width:100%;height:100%;
  border:none;pointer-events:none;
  z-index:1;
}}

/* ══════════════════════════════════════════════
   GRID
   Columns : 2vw | 36vw | 1fr | 36vw | 2vw
   Rows    : 17vh | auto | 30vh

   LEFT  panel → row 2, col 2
   RIGHT panel → row 2, col 4
   HUD visible in all other cells
   ══════════════════════════════════════════════ */
#grid{{
  position:fixed;inset:0;
  display:grid;
  grid-template-columns: 2vw 36vw 1fr 36vw 2vw;
  grid-template-rows:    17vh auto 30vh;
  pointer-events:none;   /* transparent — only panels catch events */
  z-index:500;
}}

/* ══ TITLE CELL (row 1, cols 2-4) — centred above both panels ══ */
#title-cell{{
  grid-row:1; grid-column:2/5;
  display:flex;
  align-items:flex-end;
  justify-content:center;
  padding-bottom:10px;
  pointer-events:none;
}}
#title-cell h1{{
  color:#00ffee;
  text-shadow:0 0 16px #00e5cc, 0 0 30px rgba(0,229,204,0.3);
  font-size:clamp(13px,1.4vw,20px);
  letter-spacing:5px;
  text-transform:uppercase;
  text-align:center;
  line-height:1.4;
}}
#title-cell h1 .tri{{color:#00e5cc;margin-right:8px}}

/* ══ LEFT PANEL (row 2, col 2) ══ */
#left{{
  grid-row:2; grid-column:2;
  pointer-events:all;
  background:#000000;
  border:1px solid rgba(0,229,204,0.4);
  border-radius:4px;
  padding:14px 18px 16px;
  align-self:start;
  max-height:53vh;
  overflow-y:auto;
  overflow-x:hidden;
  scrollbar-width:thin;
  scrollbar-color:#00e5cc #001813;
}}
#left::-webkit-scrollbar{{width:4px}}
#left::-webkit-scrollbar-thumb{{background:#00e5cc;border-radius:2px}}

/* ══ RIGHT PANEL (row 2, col 4) ══ */
#right{{
  grid-row:2; grid-column:4;
  pointer-events:all;
  background:#000000;
  border:1px solid rgba(0,229,204,0.4);
  border-radius:4px;
  padding:16px 18px;
  align-self:start;
  max-height:53vh;
  overflow-y:auto;
  overflow-x:hidden;
  scrollbar-width:thin;
  scrollbar-color:#00e5cc #001813;
}}
#right::-webkit-scrollbar{{width:4px}}
#right::-webkit-scrollbar-thumb{{background:#00e5cc;border-radius:2px}}

/* ══ Shared styles ══ */

/* select rows */
.sr{{display:flex;flex-direction:column;gap:3px;margin-bottom:11px}}
.sr>label{{
  font-size:9.5px;letter-spacing:1.5px;text-transform:uppercase;
  color:rgba(0,229,204,0.6);
}}
.sr select{{
  background:#000;color:#00ffee;
  border:1px solid rgba(0,229,204,0.4);border-radius:3px;
  padding:6px 26px 6px 8px;
  font-family:'Courier New',monospace;font-size:11.5px;
  width:100%;cursor:pointer;outline:none;
  -webkit-appearance:none;appearance:none;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%2300e5cc'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:right 9px center;
}}
.sr select:focus{{border-color:#00e5cc;box-shadow:0 0 7px rgba(0,229,204,0.35)}}
.sr select option{{background:#000d0d;color:#00e5cc}}

.hdiv{{border:none;border-top:1px solid rgba(0,229,204,0.18);margin:10px 0 12px}}

/* input field rows */
.frow{{margin-bottom:9px}}
.frow label{{
  display:block;font-size:10px;
  color:rgba(0,229,204,0.8);
  margin-bottom:4px;letter-spacing:0.3px;
}}
.frow input{{
  width:100%;
  background:#000;color:#00ffee;
  border:1px solid rgba(0,229,204,0.4);border-radius:3px;
  padding:6px 8px;
  font-family:'Courier New',monospace;font-size:12px;
  outline:none;
}}
.frow input:focus{{border-color:#00e5cc;box-shadow:0 0 6px rgba(0,229,204,0.35)}}
.frow input::-webkit-inner-spin-button,
.frow input::-webkit-outer-spin-button{{-webkit-appearance:none;margin:0}}
.frow input[type=number]{{-moz-appearance:textfield}}

/* caption */
.cap{{
  font-size:9px;color:rgba(0,229,204,0.4);
  letter-spacing:1px;margin:0 0 12px;
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
  transition:background .2s,box-shadow .2s;
  margin-bottom:14px;
}}
#cbtn:hover{{background:rgba(0,229,204,0.2);box-shadow:0 0 18px rgba(0,229,204,0.5)}}
#cbtn:active{{background:rgba(0,229,204,0.32)}}

/* result box — bold green */
.rbox{{
  padding:11px 14px;border-radius:3px;
  font-size:12px;line-height:1.9;
  word-break:break-word;margin-bottom:13px;
}}
.rbox.ok{{
  background:rgba(0,70,35,0.92);
  border:1px solid rgba(0,200,100,0.45);
  border-left:3px solid #00e87a;
  color:#00ffcc;
  font-weight:bold;
}}
.rbox.ok small{{
  font-weight:normal;font-size:10.5px;
  color:rgba(0,255,200,0.75);
}}
.rbox.err{{
  background:rgba(50,0,0,0.95);
  border:1px solid rgba(255,60,60,0.35);
  border-left:3px solid #ff4444;
  color:#ff9090;
}}

/* spec table */
.sp-title{{
  font-size:10px;letter-spacing:1px;
  color:#00e5cc;text-transform:uppercase;
  margin-bottom:6px;font-weight:bold;
}}
.stbl{{width:100%;border-collapse:collapse;font-size:10.5px}}
.stbl th{{
  color:#00e5cc;background:rgba(0,229,204,0.08);
  border-bottom:1px solid rgba(0,229,204,0.3);
  padding:5px 9px;text-align:left;
  letter-spacing:1px;font-size:9.5px;text-transform:uppercase;
}}
.stbl td{{
  color:rgba(0,229,204,0.85);
  border-bottom:1px solid rgba(0,229,204,0.1);
  padding:5px 9px;
}}
.stbl tr:last-child td{{border-bottom:none}}
</style>
</head>
<body>

<!-- HUD: animated full-screen background -->
{"<iframe id='hud' src='" + hud_src + "' sandbox='allow-scripts'></iframe>" if hud_src else ""}

<div id="grid">

  <!-- TITLE: centred above both panels (row1, cols 2-4) -->
  <div id="title-cell">
    <h1><span class="tri">⟁</span>HIGH-PRESSURE GAS FLOW CALCULATOR</h1>
  </div>

  <!-- LEFT PANEL: selects + fields -->
  <div id="left">
    <div class="sr">
      <label>Gas type:</label>
      <select id="gasSelect" onchange="rebuildFields()">{gas_opts}</select>
    </div>
    <div class="sr">
      <label>Calculation type:</label>
      <select id="calcSelect" onchange="rebuildFields()">{calc_opts}</select>
    </div>
    <hr class="hdiv">
    <div id="input-fields">{fields_html}</div>
  </div><!-- #left -->

  <!-- RIGHT PANEL: caption + button + result + table -->
  <div id="right">
    <div class="cap">Friction factor (f) = {FRICTION_FACTOR:.2f}&nbsp;&nbsp;·&nbsp;&nbsp;fixed constant</div>
    <button id="cbtn" onclick="doCalc()">▶ &nbsp;Calculate</button>
    <div id="result-area">{result_html}</div>
  </div>

</div><!-- #grid -->

<script>
const CF  = {CF_JS};
const DEF = {DEF_JS};

function toKey(f){{
  return f.replace(/ /g,'_').replace(/[()]/g,'').replace(/\\//g,'').replace(/°/g,'deg');
}}

function renderField(f){{
  const k = toKey(f), v = DEF[f] ?? 0;
  return `<div class="frow">
    <label>${{f}}</label>
    <input type="number" id="${{k}}" value="${{v}}" step="any">
  </div>`;
}}

function rebuildFields(){{
  const calc   = document.getElementById('calcSelect').value;
  const fields = CF[calc] || [];
  document.getElementById('input-fields').innerHTML = fields.map(renderField).join('');
  document.getElementById('result-area').innerHTML  = '';
}}

function doCalc(){{
  const gas  = document.getElementById('gasSelect').value;
  const calc = document.getElementById('calcSelect').value;
  const p    = new URLSearchParams();
  p.set('action', 'calc');
  p.set('gas',    gas);
  p.set('calc',   calc);
  (CF[calc] || []).forEach(f => {{
    const el = document.getElementById(toKey(f));
    p.set(toKey(f), el ? el.value : (DEF[f] ?? 0));
  }});
  window.parent.location.search = '?' + p.toString();
}}
</script>
</body>
</html>"""

# ── Render ────────────────────────────────────────────────────────────────────
components.html(page_html, height=800, scrolling=False)
