# Magneto-Rossby Dispersion Dashboard

An interactive Streamlit app for visualizing the sectoral dispersion curves of hydrodynamic Rossby waves, fast/slow magneto-Rossby waves, and Alfvén waves as a function of azimuthal order *m* and magnetic field strength *B*.

**Live app:** [https://magneto-rossby-dashboard.streamlit.app/](https://magneto-rossby-dashboard.streamlit.app/)

![Dashboard Example](Figures/example.png)

---

## Technical Skills

- **Web app development:** built and deployed an interactive Streamlit application accessible via public URL
- **Interactive data visualization:** real-time parameter controls with dynamic plot updates, axis zoom sliders, and multi-curve overlays, which are designed for scientific exploration rather than static presentation
- **Scientific computing:** quadratic dispersion relation solver, Alfvén speed computation, unit conversion between physical and dimensionless quantities
- **Data products:** tabular data view and one-click CSV export of all plotted curves
- **Deployment:** hosted on Streamlit Community Cloud; reproducible environment via `requirements.txt`
- **Libraries:** `streamlit`, `numpy`, `pandas`, `matplotlib`

---

## What the App Does

Users can adjust physical parameters from the sidebar and immediately see how the dispersion curves respond, no coding required. This makes it useful both as a research tool for exploring parameter space and as an educational tool for understanding wave physics.

**Adjustable parameters:**
- Rotation rate Ω, solar radius R₀, and plasma density ρ
- Azimuthal order range *m*
- One or more magnetic field strengths *B* (via presets or a custom list)
- Axis limits for zooming into any region of *m*–ν space

**Output:**
- Live dispersion curve plot with labeled wave families and color-coded *B* values
- Tabular data view of all computed curves
- CSV download of the full dataset

---

## Physics Background

The dispersion relations implemented here follow [Zaqarashvili, Oliver, Ballester & Shergelashvili 2007](https://www.aanda.org/articles/aa/abs/2007/30/aa7382-07/aa7382-07.html), with Alfvén and Hydrodynamical (HD) Rossby branches included for reference.

$$n(n+1)\,x^2 + x + \alpha^2\bigl[2 - n(n+1)\bigr] = 0, \quad x = \lambda/s$$

where:
- *s* = *n* = *m* for sectoral modes
- λ is the dimensionless eigenvalue related to mode frequency by ν = 2Ω λ
- α = v_A / (2Ω R₀) is the dimensionless magnetic parameter
- v_A = B / √(μ₀ ρ) is the Alfvén speed

The two roots correspond to the **fast** (closer to the HD Rossby branch) and **slow** magneto-Rossby branches. The **Alfvén** frequency is given by ν_A = (m / R₀) v_A / 2π. In the limit α → 0, the slow branch recovers the purely hydrodynamic Rossby dispersion relation.

Some notes: (1) density (ρ) has a large effect on the Alfvén speed (and thus α); (2) for very large B or low ρ, the slow magneto-Rossby branch can depart significantly from the HD Rossby baseline; (3) the Alfvén branches are symmetric about ν = 0; only the prograde branch is labeled in the legend to avoid duplication; and (4) all frequencies are in nHz, consistent with helioseismology conventions.


---

## Quickstart

### Run locally

```bash
pip install -r requirements.txt
streamlit run dashboard_magnetorossby.py
```

### View online

[https://magneto-rossby-dashboard.streamlit.app/](https://magneto-rossby-dashboard.streamlit.app/)

---
