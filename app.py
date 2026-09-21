import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="R-Loop Nexus", page_icon="🧬", layout="wide")

DEFAULT_THRESHOLD = 30
K_SWEEP = np.round(np.arange(0.10, 0.91, 0.10), 2)


def r_loop_probability(n, k, threshold=DEFAULT_THRESHOLD):
    n = np.asarray(n, dtype=float)
    return 1.0 / (1.0 + np.exp(-k * (n - threshold)))


def transition_width(k, p_low=0.10, p_high=0.90):
    k = np.asarray(k, dtype=float)
    return (np.log(p_high / (1.0 - p_high)) - np.log(p_low / (1.0 - p_low))) / k


def probability_table(k_values, n_points, threshold):
    data = {"Repeat Length (n)": n_points}
    for kval in k_values:
        data[f"k = {kval:.2f}"] = np.round(r_loop_probability(n_points, kval, threshold), 3)
    return pd.DataFrame(data)


st.markdown(
    """
    <div style="background:#0B1F3A;padding:18px 22px;border-radius:14px;margin-bottom:14px;">
        <div style="color:#9BD0FF;font-size:13px;letter-spacing:2px;font-weight:700;">EDUCATIONAL SIMULATION</div>
        <div style="color:white;font-size:40px;font-weight:800;line-height:1.1;">R-Loop Nexus</div>
        <div style="color:#D7E8F7;font-size:16px;margin-top:4px;">
            C9orf72 GGGGCC theoretical R-loop model
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)
st.caption("R-Loop Nexus is a teaching model. Not a diagnostic or medical test.")

st.markdown("""
This app draws one curve.

- **Left side of the curve** = short DNA repeats, low modeled chance of an R-loop
- **Middle dashed line** = 30 repeats, the center of the model
- **Right side of the curve** = long DNA repeats, high modeled chance of an R-loop

An **R-loop** is when RNA sticks to DNA instead of leaving.
R-Loop Nexus does not test a person and does not diagnose disease.
""")

c1, c2, c3 = st.columns(3)
c1.success("1. Move Repeat length")
c2.success("2. Move Steepness")
c3.success("3. Read the box below the graph")

st.sidebar.markdown("### R-Loop Nexus")
st.sidebar.caption("Controls")
if st.sidebar.button("Reset to example"):
    st.session_state["n_inspect"] = 30
    st.session_state["k"] = 0.20
    st.session_state["compare"] = True
    st.session_state["k2"] = 0.50
    st.session_state["threshold"] = 30
    st.session_state["max_n"] = 100

n_inspect = st.sidebar.slider("Repeat length to inspect (n)", 0, 100, 30, 1, key="n_inspect")
k = st.sidebar.slider("Steepness (k)", 0.10, 0.90, 0.20, 0.01, key="k")
compare = st.sidebar.checkbox("Compare a second k value", value=True, key="compare")
k2 = None
if compare:
    k2 = st.sidebar.slider("Second steepness (k2)", 0.10, 0.90, 0.50, 0.01, key="k2")
threshold = st.sidebar.slider("Center of the curve (repeats)", 20, 35, 30, 1, key="threshold")
max_n = st.sidebar.slider("Max repeat length shown", 50, 100, 100, 10, key="max_n")

st.sidebar.markdown("---")
st.sidebar.caption(
    "R-Loop Nexus example: n = 30, k = 0.20. "
    "Optional comparison k = 0.50. "
    "These are teaching values, not lab measurements."
)

n_sample = np.arange(0, max_n + 1, 10)
n_fine = np.linspace(0, max_n, 400)
P_fine = r_loop_probability(n_fine, k, threshold)
P_sample = r_loop_probability(n_sample, k, threshold)
width = float(transition_width(k))
p_now = float(r_loop_probability(n_inspect, k, threshold))
p_mid = float(r_loop_probability(threshold, k, threshold))

P2_fine = P2_sample = width2 = p2_now = None
if k2 is not None:
    P2_fine = r_loop_probability(n_fine, k2, threshold)
    P2_sample = r_loop_probability(n_sample, k2, threshold)
    width2 = float(transition_width(k2))
    p2_now = float(r_loop_probability(n_inspect, k2, threshold))

if k < 0.25:
    how = "slow and spread out"
    tip = "The curve climbs gently."
elif k < 0.55:
    how = "medium"
    tip = "The climb is easy to see around the dashed line."
else:
    how = "sudden"
    tip = "The curve jumps up in a short span."

if n_inspect < threshold - 5:
    zone = "You are on the left side of the curve, below the center."
    zone_meaning = "At this shorter repeat length, the modeled R-loop likelihood is still low."
elif abs(n_inspect - threshold) <= 5:
    zone = "You are near the center of the curve."
    zone_meaning = (
        "The model is built so the middle is 50% at this repeat length. "
        "That 50% is a modeling choice, not a lab-measured chance of disease."
    )
else:
    zone = "You are on the right side of the curve, above the center."
    zone_meaning = "At this longer repeat length, the modeled R-loop likelihood is high."

col1, col2 = st.columns([1.35, 1], gap="large")

with col1:
    fig, ax = plt.subplots(figsize=(8, 5.4))
    ax.plot(n_fine, P_fine, label=f"k = {k:.2f}", color="steelblue", linewidth=2.4)
    ax.scatter(n_sample, P_sample, color="steelblue", zorder=5)
    if k2 is not None:
        ax.plot(n_fine, P2_fine, label=f"k = {k2:.2f}", color="firebrick", linewidth=2.4)
        ax.scatter(n_sample, P2_sample, color="firebrick", zorder=5)
    ax.axvline(threshold, color="gray", linestyle="--", linewidth=1.2, label=f"Center (n = {threshold})")
    ax.axhline(0.10, color="lightgray", linestyle=":", linewidth=1)
    ax.axhline(0.50, color="lightgray", linestyle=":", linewidth=1)
    ax.axhline(0.90, color="lightgray", linestyle=":", linewidth=1)
    ax.axvline(n_inspect, color="#2E8B57", linestyle="-.", linewidth=1.2, label=f"Inspect n = {n_inspect}")
    ax.scatter([n_inspect], [p_now], color="#2E8B57", s=70, zorder=6)
    ax.set_xlabel("GGGGCC repeat length (n)")
    ax.set_ylabel("Model output (0 to 1)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("R-Loop Nexus")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.28)
    st.pyplot(fig)
    st.caption("Green marker = the repeat length you are inspecting. Dashed line = center of the model.")

with col2:
    st.subheader("What this means right now")
    st.metric("Repeat length", f"{n_inspect}")
    st.metric("Model output", f"{p_now:.3f}")
    st.metric("10% to 90% width", f"{width:.1f} repeats")
    st.write(zone)
    st.write(zone_meaning)

st.info(f"""
**R-Loop Nexus settings right now**

- Center of the curve: **{threshold} repeats**
- Steepness k: **{k:.2f}** ({how})
- Value at the center: **{p_mid:.0%}**
- Distance from 10% to 90% on the curve: **{width:.1f} repeats**

{tip}

Try this: raise Steepness and watch the curve get taller faster near the dashed line.
""")

if k2 is not None:
    st.write(
        f"Comparison: at n = {n_inspect}, the first curve is {p_now:.3f} and the second curve "
        f"(k = {k2:.2f}) is {p2_now:.3f}. The 10% to 90% width changes from {width:.1f} to {width2:.1f} repeats."
    )

st.markdown("### Why R-Loop Nexus is useful")
u1, u2, u3 = st.columns(3)
with u1:
    st.markdown("**Students**")
    st.write("See how a threshold curve works: low, then a rise, then high.")
with u2:
    st.markdown("**Teachers**")
    st.write("Show a genetics idea without laboratory equipment.")
with u3:
    st.markdown("**Everyone else**")
    st.write("Explore the graph. Do not use it as a medical test.")

tab1, tab2, tab3, tab4 = st.tabs(
    [
        "Simple terms",
        "P(n) for selected k",
        "Table 2 — k-sweep",
        "Table 3 — transition width",
    ]
)

with tab1:
    st.markdown("""
**n (repeat length)**  
How many times GGGGCC is copied. Bigger n = farther right on the graph.

**k (steepness)**  
How quickly the curve rises. Bigger k = sharper rise.

**P(n)**  
The height of the curve. 0 is low, 1 is high. This is a model output, not a lab result.

**Why 30?**  
Researchers often use about 30 C9orf72 repeats as a reference point.
R-Loop Nexus places the middle of the curve there. That does **not** mean a real R-loop is exactly 50% likely at 30 repeats.

**Transition width**  
The distance between 10% and 90% on the curve. In this formula it equals 2 × ln(9) / k.

**What R-Loop Nexus cannot do**  
It cannot diagnose ALS or FTD, read a person's DNA, or replace genetic testing.
    """)

with tab2:
    selected_ks = [k] if k2 is None else [k, k2]
    table_selected = probability_table(selected_ks, n_sample, threshold)
    st.write("Model outputs at every 10 repeats for the k value(s) you selected.")
    st.dataframe(table_selected, use_container_width=True, hide_index=True)
    st.download_button(
        "Download selected-k table (CSV)",
        table_selected.to_csv(index=False).encode("utf-8"),
        "rloop_nexus_selected_k.csv",
        "text/csv",
    )

with tab3:
    st.write("Predicted P(n) at 10-repeat intervals for k = 0.10 to 0.90. Every column is 0.500 at the center.")
    table_sweep = probability_table(K_SWEEP, n_sample, threshold)
    st.dataframe(table_sweep, use_container_width=True, hide_index=True)
    st.download_button(
        "Download k-sweep table (CSV)",
        table_sweep.to_csv(index=False).encode("utf-8"),
        "rloop_nexus_ksweep.csv",
        "text/csv",
    )

with tab4:
    widths = np.round(transition_width(K_SWEEP), 2)
    width_row = {"k": "Transition width (n at P=0.90 minus n at P=0.10)"}
    for k_i, w in zip(K_SWEEP, widths):
        width_row[f"{k_i:.2f}"] = w
    table_width = pd.DataFrame([width_row])
    st.write("Larger k produces a narrower transition. This is a property of the formula, not a lab measurement.")
    st.dataframe(table_width, use_container_width=True, hide_index=True)
    st.download_button(
        "Download transition-width table (CSV)",
        table_width.to_csv(index=False).encode("utf-8"),
        "rloop_nexus_transition_width.csv",
        "text/csv",
    )

st.markdown("---")
st.markdown("**R-Loop Nexus** · Educational computational model · Not for diagnosis")
st.caption("Parameters are teaching assumptions, not experimentally calibrated constants.")
