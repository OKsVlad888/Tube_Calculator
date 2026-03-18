import math
import base64
from pathlib import Path

import streamlit as st

# Page configuration (kept as in your code)
st.set_page_config(
    page_title="High-Pressure Gas Flow Calculator",
    page_icon="💨",
    layout="centered",
)

# ── Load the HUD background HTML (UNCHANGED BEHAVIOR) ────────────────────────

def load_hud_html() -> str | None:
    path = Path(__file__).resolve().parent / "hud_background.html"
    if not path.exists():
        return None
    return path.read_text(encoding="utf-8")


hud_html_raw = load_hud_html()


def make_iframe_src(html: str) -> str:
    encoded = base64.b64encode(html.encode("utf-8")).decode("utf-8")
    return f"data:text/html;base64,{encoded}"


# ── Inject CSS that embeds the HUD iframe as a fixed full-screen background ──
# NOTE: We keep the HUD layer logic intact. We only add a scrollable inner panel
#       for inputs+outputs so the page itself does NOT scroll.

def inject_styles(iframe_src: str | None) -> None:
    bg_iframe = (
        f"""
        <iframe
            src=\"{iframe_src}\"
            style=\"position: fixed; top: 0; left: 0; width: 100vw; height: 100vh; border: none; z-index: 0;\"
        ></iframe>
        """
        if iframe_src
        else ""
    )

    css = """
    <style>
      /* Make the Streamlit app transparent and keep content above HUD */
      .stApp { background: transparent; }
      .stApp > header { background: transparent; }

      /* Ensure all Streamlit content is above the HUD iframe */
      .stApp > div { position: relative; z-index: 2; }

      /* --- NEW: Create a fixed-height internal scroll panel --- */
      /* Centered layout: constrain width similar to your current panel */
      .calc-panel {
        background: rgba(0, 0, 0, 0.62);
        border: 1px solid rgba(0, 255, 240, 0.20);
        border-radius: 14px;
        padding: 16px 18px;
        box-shadow: 0 0 0 1px rgba(0,255,240,0.06), 0 18px 40px rgba(0,0,0,0.55);
        backdrop-filter: blur(6px);
      }

      /* Reduce extra whitespace so the panel feels compact */
      .calc-panel .stMarkdown, .calc-panel label, .calc-panel p { color: #bffdf7; }

      /* IMPORTANT: prevent overall page scrolling; scroll should be inside container */
      html, body { overflow: hidden; }

      /* Streamlit main area can still scroll if something forces it; limit it */
      section.main { overflow: hidden; }
    </style>
    """

    st.markdown(bg_iframe + css, unsafe_allow_html=True)


iframe_src = make_iframe_src(hud_html_raw) if hud_html_raw else None
inject_styles(iframe_src)

if not hud_html_raw:
    st.warning("⚠️ 'hud_background.html' not found – using plain dark background.")

# ═════════════════════════════════════════════════════════════════════════════
# GAS DATA & CONSTANTS (UNCHANGED)
# ═════════════════════════════════════════════════════════════════════════════

GAS_DATA = {
    "N2": 0.028013,
    "O2": 0.031999,
    "Ar": 0.039948,
    "CO2": 0.04401,
    "He": 0.0040026,
    "H2": 0.002016,
    "CH4": 0.01604,
    "C2H2": 0.02604,
    "Forming Gas 1": 0.03881,
    "Forming Gas 2": 0.02671,
    "Air": 0.02897,
}

GAS_DISPLAY_NAMES = {
    "N2": "N2 (Nitrogen)",
    "O2": "O2 (Oxygen)",
    "Ar": "Ar (Argon)",
    "CO2": "CO2 (Carbon Dioxide)",
    "He": "He (Helium)",
    "H2": "H2 (Hydrogen)",
    "CH4": "CH4 (Methane)",
    "C2H2": "C2H2 (Acetylene)",
    "Forming Gas 1": "(H2-3% + Ar-97%) Forming Gas 1",
    "Forming Gas 2": "(H2-5% + N2-95%) Forming Gas 2",
    "Air": "Air (Dry Air)",
}

CALC_FIELDS = {
    "Pipe Diameter (mm)": [
        "Temperature (°C)",
        "Inlet Pressure (bar)",
        "Outlet Pressure (bar)",
        "Pipe Length (m)",
        "Flow Rate (LPM)",
    ],
    "Flow Rate (LPM)": [
        "Temperature (°C)",
        "Inlet Pressure (bar)",
        "Outlet Pressure (bar)",
        "Pipe Length (m)",
        "Pipe Diameter (mm)",
    ],
    "Pipe Length (m)": [
        "Temperature (°C)",
        "Inlet Pressure (bar)",
        "Outlet Pressure (bar)",
        "Pipe Diameter (mm)",
        "Flow Rate (LPM)",
    ],
    "Inlet Pressure (bar)": [
        "Temperature (°C)",
        "Outlet Pressure (bar)",
        "Pipe Length (m)",
        "Pipe Diameter (mm)",
        "Flow Rate (LPM)",
    ],
    "Outlet Pressure (bar)": [
        "Temperature (°C)",
        "Inlet Pressure (bar)",
        "Pipe Length (m)",
        "Pipe Diameter (mm)",
        "Flow Rate (LPM)",
    ],
}

FRICTION_FACTOR = 0.02

# ═════════════════════════════════════════════════════════════════════════════
# CALCULATION FUNCTIONS (UNCHANGED)
# ═════════════════════════════════════════════════════════════════════════════


def ideal_gas_density(P_pa: float, T_k: float, gas: str) -> float:
    R = 8.314
    M = GAS_DATA[gas]
    return (P_pa * M) / (R * T_k)


def calc_required_diameter(P_in, P_out, T_c, L, Q, gas):
    P_in_pa = P_in * 1e5
    P_out_pa = P_out * 1e5
    T_k = T_c + 273.15
    Q_m3_s = Q / 1000.0 / 60.0
    rho_avg = (
        ideal_gas_density(P_in_pa, T_k, gas)
        + ideal_gas_density(P_out_pa, T_k, gas)
    ) / 2.0
    delta_p = P_in_pa - P_out_pa
    if delta_p <= 0:
        raise ValueError("Inlet pressure must be greater than outlet pressure.")
    D_m = (
        (FRICTION_FACTOR * L * 8.0 * rho_avg * Q_m3_s**2)
        / (math.pi**2 * delta_p)
    ) ** 0.2
    return D_m * 1000.0


def calc_flow_rate(P_in, P_out, T_c, L, D, gas):
    P_in_pa = P_in * 1e5
    P_out_pa = P_out * 1e5
    T_k = T_c + 273.15
    D_m = D / 1000.0
    rho_avg = (
        ideal_gas_density(P_in_pa, T_k, gas)
        + ideal_gas_density(P_out_pa, T_k, gas)
    ) / 2.0
    delta_p = P_in_pa - P_out_pa
    if delta_p <= 0:
        raise ValueError("Inlet pressure must be higher than outlet pressure.")
    Q_m3_s = math.sqrt(
        (delta_p * math.pi**2 * D_m**5) / (8.0 * FRICTION_FACTOR * L * rho_avg)
    )
    return Q_m3_s * 1000.0 * 60.0


def calc_max_length(P_in, P_out, T_c, D, Q, gas):
    P_in_pa = P_in * 1e5
    P_out_pa = P_out * 1e5
    T_k = T_c + 273.15
    D_m = D / 1000.0
    Q_m3_s = Q / 1000.0 / 60.0
    rho_avg = (
        ideal_gas_density(P_in_pa, T_k, gas)
        + ideal_gas_density(P_out_pa, T_k, gas)
    ) / 2.0
    delta_p = P_in_pa - P_out_pa
    if delta_p <= 0:
        raise ValueError("Inlet pressure must be greater than outlet pressure.")
    return (
        delta_p * math.pi**2 * D_m**5
    ) / (8.0 * FRICTION_FACTOR * rho_avg * Q_m3_s**2)


def calc_outlet_pressure(P_in, T_c, L, D, Q, gas):
    P_in_pa = P_in * 1e5
    T_k = T_c + 273.15
    D_m = D / 1000.0
    Q_m3_s = Q / 1000.0 / 60.0
    P_out_guess = P_in
    for _ in range(20):
        P_out_pa = P_out_guess * 1e5
        rho_avg = (
            ideal_gas_density(P_in_pa, T_k, gas)
            + ideal_gas_density(P_out_pa, T_k, gas)
        ) / 2.0
        delta_p = (
            8.0 * FRICTION_FACTOR * L * rho_avg * Q_m3_s**2
        ) / (math.pi**2 * D_m**5)
        P_out_calc = P_in - delta_p / 1e5
        if abs(P_out_calc - P_out_guess) < 0.001:
            return max(P_out_calc, 0.0)
        P_out_guess = P_out_calc
    return max(P_out_guess, 0.0)


def calc_required_inlet_pressure(P_out, T_c, L, D, Q, gas):
    P_low, P_high = P_out, P_out + 100.0
    while P_high < P_out + 2000:
        if calc_outlet_pressure(P_high, T_c, L, D, Q, gas) <= P_out:
            break
        P_high += 100.0
    for _ in range(50):
        P_mid = (P_low + P_high) / 2.0
        P_out_mid = calc_outlet_pressure(P_mid, T_c, L, D, Q, gas)
        if abs(P_out_mid - P_out) < 0.01:
            return P_mid
        if P_out_mid > P_out:
            P_low = P_mid
        else:
            P_high = P_mid
    return (P_low + P_high) / 2.0


def determine_pipe_options(D, P, gas):
    if gas == "O2":
        return ["1\" tube (Spec S14) – required for O2"], "1\" tube (Spec S14)"

    opts, rec = [], None
    if D <= 4.0:
        if P <= 200:
            opts = [
                "1/4\" tube (Spec S6) – up to 200 bar",
                "1/4\" tube (Spec S9) – up to 1379 bar",
                "1/4\" tube (Spec S12) – above 1379 bar",
            ]
            rec = "1/4\" tube (Spec S6)"
        elif P <= 1379:
            opts = [
                "1/4\" tube (Spec S9) – up to 1379 bar",
                "1/4\" tube (Spec S12) – above 1379 bar",
            ]
            rec = "1/4\" tube (Spec S9)"
        else:
            opts = ["1/4\" tube (Spec S12) – above 1379 bar"]
            rec = "1/4\" tube (Spec S12)"
    elif D <= 7.0:
        if P <= 140:
            opts = [
                "3/8\" tube (Spec S16) – up to 140 bar",
                "3/8\" tube (Spec S9) – above 140 bar",
            ]
            rec = "3/8\" tube (Spec S16)"
        else:
            opts = ["3/8\" tube (Spec S9) – above 140 bar"]
            rec = "3/8\" tube (Spec S9)"
    elif D <= 21.0:
        if P <= 20:
            opts = [
                "3/4\" tube (Spec S15) – up to 20 bar",
                "1\" tube (Spec S14) – above 20 bar",
            ]
            rec = "3/4\" tube (Spec S15)"
        else:
            opts = ["1\" tube (Spec S14) – above 20 bar"]
            rec = "1\" tube (Spec S14)"
    else:
        opts = ["Special piping required (outside standard range)"]
        rec = "Special piping (outside range)"

    return opts, rec


# ═════════════════════════════════════════════════════════════════════════════
# UI (UPDATED: internal scroll container for the whole calculator panel)
# ═════════════════════════════════════════════════════════════════════════════

# We render the panel wrapper, then put all UI inside a fixed-height container.
st.markdown('<div class="calc-panel">', unsafe_allow_html=True)

# Choose a height that fits most screens; users will scroll INSIDE this.
# You can tune 700/750 based on your preference.
with st.container(height=740):
    st.title("⟁ HIGH-PRESSURE GAS FLOW CALCULATOR")

    gas_type = st.selectbox(
        "Gas type:",
        list(GAS_DATA.keys()),
        format_func=lambda x: GAS_DISPLAY_NAMES[x],
    )
    calc_type = st.selectbox("Calculation type:", list(CALC_FIELDS.keys()))

    DEFAULTS = {
        "Temperature (°C)": 25.0,
        "Inlet Pressure (bar)": 100.0,
        "Outlet Pressure (bar)": 10.0,
        "Pipe Length (m)": 10.0,
        "Pipe Diameter (mm)": 10.0,
        "Flow Rate (LPM)": 100.0,
    }

    values: dict[str, float] = {}
    for field in CALC_FIELDS[calc_type]:
        values[field] = st.number_input(field, value=float(DEFAULTS.get(field, 0.0)))

    st.caption(f"Friction factor (f) = {FRICTION_FACTOR:.2f} • fixed constant")

    # Session state
    if "result" not in st.session_state:
        st.session_state.update(result=None, error=False, specs=[])

    if st.button("▶ CALCULATE"):
        try:
            T_c = float(values.get("Temperature (°C)", 0.0))
            P_in = float(values.get("Inlet Pressure (bar)", 0.0))
            P_out = float(values.get("Outlet Pressure (bar)", 0.0))
            L_val = float(values.get("Pipe Length (m)", 0.0))
            D_val = float(values.get("Pipe Diameter (mm)", 0.0))
            Q_val = float(values.get("Flow Rate (LPM)", 0.0))

            msg = ""
            options: list[str] = []

            if calc_type == "Pipe Diameter (mm)":
                D_req = calc_required_diameter(P_in, P_out, T_c, L_val, Q_val, gas_type)
                options, rec = determine_pipe_options(D_req, max(P_in, P_out), gas_type)
                msg = f"Required Diameter: **{D_req:.2f} mm**\nRecommended: **{rec}**"

            elif calc_type == "Flow Rate (LPM)":
                Q_max = calc_flow_rate(P_in, P_out, T_c, L_val, D_val, gas_type)
                options, rec = determine_pipe_options(D_val, max(P_in, P_out), gas_type)
                msg = f"Maximum Flow Rate: **{Q_max:.1f} L/min**\nRecommended: **{rec}**"

            elif calc_type == "Pipe Length (m)":
                L_max = calc_max_length(P_in, P_out, T_c, D_val, Q_val, gas_type)
                options, rec = determine_pipe_options(D_val, max(P_in, P_out), gas_type)
                msg = f"Maximum Pipe Length: **{L_max:.1f} m**\nRecommended: **{rec}**"

            elif calc_type == "Inlet Pressure (bar)":
                P_req = calc_required_inlet_pressure(P_out, T_c, L_val, D_val, Q_val, gas_type)
                options, rec = determine_pipe_options(D_val, max(P_req, P_out), gas_type)
                msg = f"Required Inlet Pressure: **{P_req:.2f} bar**\nRecommended: **{rec}**"

            elif calc_type == "Outlet Pressure (bar)":
                P_est = calc_outlet_pressure(P_in, T_c, L_val, D_val, Q_val, gas_type)
                options, rec = determine_pipe_options(D_val, P_in, gas_type)
                msg = f"Estimated Outlet Pressure: **{P_est:.2f} bar**\nRecommended: **{rec}**"

            else:
                msg = "Error: Unsupported calculation type."

            st.session_state.update(result=msg, error=False, specs=options)

        except Exception as e:
            st.session_state.update(result=f"Error: {e}", error=True, specs=[])

    # Output (still inside the scroll panel)
    if st.session_state.result is not None:
        if st.session_state.error:
            st.error(st.session_state.result)
        else:
            st.success(st.session_state.result)

        if st.session_state.specs:
            st.markdown("**Possible pipe specifications:**")
            data = []
            for spec in st.session_state.specs:
                if "–" in spec:
                    name, detail = spec.split("–", 1)
                    data.append({"Pipe Spec": name.strip(), "Details": detail.strip()})
                else:
                    data.append({"Pipe Spec": spec, "Details": ""})
            st.table(data)

st.markdown('</div>', unsafe_allow_html=True)
