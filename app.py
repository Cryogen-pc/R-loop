import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="R-loop Formation Probability Model", layout="wide")

# ---------------------------------------------------------
# Core model
# ---------------------------------------------------------
def r_loop_probability(n, k, threshold=30):
    return 1 / (1 + np.exp(-k * (n - threshold)))

st.title("R-loop Formation Probability Model")
st.caption(
    "A computational logistic model of R-loop formation probability in C9orf72 "
    "GGGGCC repeat expansions."
)

# ---------------------------------------------------------
# Sidebar controls
# ---------------------------------------------------------
st.sidebar.header("Model Parameters")

k_base = st.sidebar.slider(
    "Baseline steepness (k_base)", min_value=0.05, max_value=1.0, value=0.20, step=0.01
)
k_therm = st.sidebar.slider(
    "Thermodynamic-adjusted steepness (k_therm)", min_value=0.05, max_value=1.0, value=0.50, step=0.01
)
threshold = st.sidebar.slider(
    "Pathogenicity threshold (repeats)", min_value=10, max_value=60, value=30, step=1
)
max_n = st.sidebar.slider(
    "Max repeat length shown (n)", min_value=50, max_value=300, value=100, step=10
)

st.sidebar.markdown("---")
st.sidebar.caption(
    "k_base = 0.2 and k_therm = 0.5 with threshold = 30 reproduce the study's "
    "original baseline and thermodynamic-adjusted scenarios."
)

# ---------------------------------------------------------
# Data
# ---------------------------------------------------------
n_sample = np.arange(0, max_n + 1, 10)
P_base_sample = r_loop_probability(n_sample, k_base, threshold)
P_therm_sample = r_loop_probability(n_sample, k_therm, threshold)

n_fine = np.linspace(0, max_n, 400)
P_base_fine = r_loop_probability(n_fine, k_base, threshold)
P_therm_fine = r_loop_probability(n_fine, k_therm, threshold)

pct_diff = np.where(
    P_base_sample != 0,
    (P_therm_sample - P_base_sample) / P_base_sample * 100,
    np.nan,
)

# ---------------------------------------------------------
# Layout: plot + stats
# ---------------------------------------------------------
col1, col2 = st.columns([2, 1])

with col1:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(n_fine, P_base_fine, label=f"Baseline (k={k_base:.2f})", color="steelblue", linewidth=2)
    ax.plot(n_fine, P_therm_fine, label=f"Thermodynamic-adjusted (k={k_therm:.2f})", color="firebrick", linewidth=2)
    ax.axvline(x=threshold, color="gray", linestyle="--", linewidth=1, label=f"Threshold (n={threshold})")
    ax.scatter(n_sample, P_base_sample, color="steelblue", zorder=5)
    ax.scatter(n_sample, P_therm_sample, color="firebrick", zorder=5)
    ax.set_xlabel("GGGGCC Repeat Length (n)")
    ax.set_ylabel("Predicted R-loop Formation Probability, P(n)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("R-loop Formation Probability vs. Repeat Length")
    ax.legend()
    ax.grid(alpha=0.3)
    st.pyplot(fig)

with col2:
    st.subheader("Summary Statistics")
    st.metric("Mean P(n) — Baseline", f"{np.mean(P_base_sample):.3f}")
    st.metric("Mean P(n) — Thermo-adjusted", f"{np.mean(P_therm_sample):.3f}")
    st.metric("Std Dev — Baseline", f"{np.std(P_base_sample):.3f}")
    st.metric("Std Dev — Thermo-adjusted", f"{np.std(P_therm_sample):.3f}")
    st.metric("P(threshold)", f"{r_loop_probability(threshold, k_base, threshold):.3f}")

# ---------------------------------------------------------
# Results table
# ---------------------------------------------------------
st.subheader("Results Table")
results_table = pd.DataFrame({
    "Repeat Length (n)": n_sample,
    "P(n) — Baseline": np.round(P_base_sample, 4),
    "P(n) — Thermo-adjusted": np.round(P_therm_sample, 4),
    "% Difference": np.round(pct_diff, 2),
})
st.dataframe(results_table, use_container_width=True, hide_index=True)

csv = results_table.to_csv(index=False).encode("utf-8")
st.download_button("Download results as CSV", csv, "rloop_results_table.csv", "text/csv")

st.markdown("---")
st.caption(
    "This is an illustrative, non-diagnostic computational model. Parameters are based on "
    "literature-informed assumptions, not experimentally calibrated constants, and should not "
    "be used for clinical prediction."
)
