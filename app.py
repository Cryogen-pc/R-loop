import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="R-loop Formation Probability Model", layout="wide")

# ---------------------------------------------------------
# Core model
# ---------------------------------------------------------
DEFAULT_THRESHOLD = 30
K_SWEEP = np.round(np.arange(0.10, 0.91, 0.10), 2)


def r_loop_probability(n, k, threshold=DEFAULT_THRESHOLD):
    n = np.asarray(n, dtype=float)
    return 1.0 / (1.0 + np.exp(-k * (n - threshold)))


def transition_width(k, p_low=0.10, p_high=0.90):
    """Repeat-length distance between P=0.10 and P=0.90."""
    k = np.asarray(k, dtype=float)
    return (np.log(p_high / (1.0 - p_high)) - np.log(p_low / (1.0 - p_low))) / k


def probability_table(k_values, n_points, threshold):
    data = {"Repeat Length (n)": n_points}
    for k in k_values:
        data[f"k = {k:.2f}"] = np.round(r_loop_probability(n_points, k, threshold), 3)
    return pd.DataFrame(data)


st.title("R-loop Formation Probability Model")
st.caption(
    "A computational logistic model of R-loop formation probability in C9orf72 "
    "GGGGCC repeat expansions. Educational / exploratory only — not a diagnostic tool."
)

# ---------------------------------------------------------
# Sidebar controls
# ---------------------------------------------------------
st.sidebar.header("Model Parameters")

k = st.sidebar.slider(
    "Steepness (k)",
    min_value=0.10,
    max_value=0.90,
    value=0.20,
    step=0.01,
    help="Larger k makes a sharper transition around the threshold.",
)

compare = st.sidebar.checkbox("Compare a second k value", value=False)
k2 = None
if compare:
    k2 = st.sidebar.slider(
        "Second steepness (k₂)",
        min_value=0.10,
        max_value=0.90,
        value=0.50,
        step=0.01,
    )

threshold = st.sidebar.slider(
    "Pathogenicity threshold (repeats)",
    min_value=20,
    max_value=35,
    value=30,
    step=1,
)

max_n = st.sidebar.slider(
    "Max repeat length shown (n)",
    min_value=50,
    max_value=100,
    value=100,
    step=10,
)

st.sidebar.markdown("---")
st.sidebar.caption(
    "The study's main analysis is a sweep of k from 0.10 to 0.90. "
    "k = 0.20 is the worked example used in Table 1. "
    "A second k is optional for comparison only; it is not a separate "
    "calibrated thermodynamic model."
)

# ---------------------------------------------------------
# Data
# ---------------------------------------------------------
n_sample = np.arange(0, max_n + 1, 10)
n_fine = np.linspace(0, max_n, 400)

P_fine = r_loop_probability(n_fine, k, threshold)
P_sample = r_loop_probability(n_sample, k, threshold)
width = float(transition_width(k))

P2_fine = P2_sample = width2 = None
if k2 is not None:
    P2_fine = r_loop_probability(n_fine, k2, threshold)
    P2_sample = r_loop_probability(n_sample, k2, threshold)
    width2 = float(transition_width(k2))

# ---------------------------------------------------------
# Layout: plot + stats
# ---------------------------------------------------------
col1, col2 = st.columns([2, 1])

with col1:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(n_fine, P_fine, label=f"k = {k:.2f}", color="steelblue", linewidth=2)
    ax.scatter(n_sample, P_sample, color="steelblue", zorder=5)
    if k2 is not None:
        ax.plot(n_fine, P2_fine, label=f"k = {k2:.2f}", color="firebrick", linewidth=2)
        ax.scatter(n_sample, P2_sample, color="firebrick", zorder=5)
    ax.axvline(
        x=threshold,
        color="gray",
        linestyle="--",
        linewidth=1,
        label=f"Threshold (n = {threshold})",
    )
    ax.axhline(0.10, color="lightgray", linestyle=":", linewidth=1)
    ax.axhline(0.90, color="lightgray", linestyle=":", linewidth=1)
    ax.set_xlabel("GGGGCC Repeat Length (n)")
    ax.set_ylabel("Predicted R-loop Formation Probability, P(n)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("R-loop Formation Probability vs. Repeat Length")
    ax.legend()
    ax.grid(alpha=0.3)
    st.pyplot(fig)

with col2:
    st.subheader("Selected k")
    st.metric("Steepness k", f"{k:.2f}")
    st.metric("P(threshold)", f"{float(r_loop_probability(threshold, k, threshold)):.3f}")
    st.metric("Transition width", f"{width:.2f} repeats")
    st.caption("Width = n at P = 0.90 minus n at P = 0.10.")
    if k2 is not None:
        st.markdown("---")
        st.metric("Second k", f"{k2:.2f}")
        st.metric("Second transition width", f"{width2:.2f} repeats")
        st.metric("Width change (k₂ − k)", f"{width2 - width:+.2f}")

# ---------------------------------------------------------
# Tabs: selected k, full sweep (Table 2), transition width (Table 3)
# ---------------------------------------------------------
tab1, tab2, tab3 = st.tabs(
    [
        "P(n) for selected k",
        "Table 2 — k-sweep of P(n)",
        "Table 3 — transition width",
    ]
)

with tab1:
    selected_ks = [k] if k2 is None else [k, k2]
    table_selected = probability_table(selected_ks, n_sample, threshold)
    st.dataframe(table_selected, use_container_width=True, hide_index=True)
    st.download_button(
        "Download selected-k table (CSV)",
        table_selected.to_csv(index=False).encode("utf-8"),
        "rloop_selected_k.csv",
        "text/csv",
    )

with tab2:
    st.write(
        "Predicted P(n) at 10-repeat intervals for the study sweep "
        "k = 0.10 to 0.90. Every curve equals 0.500 at the threshold."
    )
    table_sweep = probability_table(K_SWEEP, n_sample, threshold)
    st.dataframe(table_sweep, use_container_width=True, hide_index=True)
    st.download_button(
        "Download k-sweep table (CSV)",
        table_sweep.to_csv(index=False).encode("utf-8"),
        "rloop_ksweep_table.csv",
        "text/csv",
    )

with tab3:
    widths = np.round(transition_width(K_SWEEP), 2)
    width_row = {"k": "Transition width (n at P=0.90 − n at P=0.10)"}
    for k_i, w in zip(K_SWEEP, widths):
        width_row[f"{k_i:.2f}"] = w
    table_width = pd.DataFrame([width_row])
    st.write(
        "Larger k produces a narrower transition around the threshold. "
        "This follows the logistic relationship width ∝ 1/k, so it describes "
        "the model rather than a laboratory measurement."
    )
    st.dataframe(table_width, use_container_width=True, hide_index=True)
    st.download_button(
        "Download transition-width table (CSV)",
        table_width.to_csv(index=False).encode("utf-8"),
        "rloop_transition_width.csv",
        "text/csv",
    )

st.markdown("---")
st.caption(
    "This is an illustrative, non-diagnostic computational model. Parameters are based on "
    "literature-informed assumptions, not experimentally calibrated constants, and should not "
    "be used for clinical prediction."
)
