# 🪳 Discrete LA Model for *Tribolium* Beetle Population

This repository presents an analysis of the discrete **Larvae–Adult (LA)** population model for *Tribolium* flour beetles. The code and data support the manuscript titled:

**"Global Dynamics of a Discrete Two-Population Model for Flour Beetle Growth"**

---

## 📁 Data

The experimental data used for model fitting is provided in:  TriboliumData_2.xlsx 

---

## 🧮 MATLAB Code Overview

### `LA_fit.m`
- **Purpose:** Fits the LA model to experimental data.
- **Output:** Estimates of optimal parameters for the LA model.

### `model_simulation.m`
- **Purpose:** Simulates the LA model using the best-fit parameters from `LA_fit.m`.
- **Features:** Demonstrates convergence to equilibrium.

### `delay_eqn_b.m` and `delay_eqn_aa.m`
- **Purpose:** Solve the LA model with varied conditions for:
  - `b` (larval population)
  - `A(0)` (initial adult population)
- **Parameters:** Uses the best-fit parameters from `LA_fit.m`.

---

## 🛠 Requirements

- MATLAB (R2022a or later recommended)
- Statistics and Machine learning toolbox, Global Optimization Toolbox, Optimization Toolbox.


