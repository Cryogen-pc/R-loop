# R-loop

A Python simulation of R-loop formation probability in C9orf72 GGGGCC repeat expansions, using a logistic (sigmoid) function anchored to the clinically recognized 30-repeat pathogenicity threshold.

## Overview

This project models the probability of R-loop formation, P(n), as a function of repeat length (n), comparing two scenarios:

- **Baseline model** (k = 0.2) — standard steepness
- **Thermodynamic-adjusted model** (k = 0.5) — accounts for the increased hybrid stability of G-C-rich RNA:DNA hybrids

The model is illustrative and educational — it is **not** a diagnostic or clinical prediction tool.

## Live demo

Try the interactive app here: **[link once deployed on Streamlit Cloud]**

## Running locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Files

- `app.py` — Streamlit app with interactive sliders for steepness (k), threshold, and repeat length range
- `requirements.txt` — Python dependencies

## Model

P(n) = 1 / (1 + e^(−k(n − threshold)))

Where:
- `n` = GGGGCC repeat length
- `k` = steepness parameter
- `threshold` = pathogenicity cutoff (default: 30 repeats)

## Background

This tool was developed as part of an undergraduate research study on computational modeling of R-loop formation in C9orf72-related ALS/FTD, motivated by limited access to genetic testing infrastructure in resource-limited settings.

## Disclaimer

Parameters are based on literature-informed assumptions, not experimentally calibrated constants. This tool is intended for educational and exploratory purposes only and should not be used for clinical diagnosis or decision-making.
