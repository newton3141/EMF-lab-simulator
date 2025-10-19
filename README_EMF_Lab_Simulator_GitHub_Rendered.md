# 🧲 EMF Lab Simulator (v4.1, GitHub Rendered)
**Electronic Induction and Lenz Force Numerical Simulator**  
*Developer: Aiden Jung 정윤서 (Seoul Science High School, 2025)*

---

## 1. Overview
Numerical simulation of **Faraday’s induction experiment**: a neodymium magnet falling through a solenoid coil.
Compares induced EMF and Lenz force between theory and experiment; ideal for performance evaluation reports.

---

## 2. Core Physics Equations

### (1) Faraday’s Law
![formula](https://latex.codecogs.com/svg.latex?%5Ctext%7BEMF%7D%20%3D%20-N%20%5Cfrac%7Bd%5CPhi%7D%7Bdt%7D)
![formula](https://latex.codecogs.com/svg.latex?%5CPhi%20%3D%20%5Cint%20B_z%20%5C%2C%20dA%20%5Capprox%20A%20B_z%28z%29%20%5CRightarrow%20%5Ctext%7BEMF%7D%20%3D%20-N%20A%20%5Cfrac%7BdB_z%7D%7Bdz%7D%20v)

### (2) Finite Cylinder Magnet Field (on-axis)
![formula](https://latex.codecogs.com/svg.latex?B_z%28z%29%20%3D%20%5Cfrac%7BB_r%7D%7B2%7D%5Cleft%28%5Cfrac%7Bz%2B%5Cfrac%7BL%7D%7B2%7D%7D%7B%5Csqrt%7B%28z%2B%5Cfrac%7BL%7D%7B2%7D%29%5E2%2BR%5E2%7D%7D%20-%20%5Cfrac%7Bz-%5Cfrac%7BL%7D%7B2%7D%7D%7B%5Csqrt%7B%28z-%5Cfrac%7BL%7D%7B2%7D%29%5E2%2BR%5E2%7D%7D%5Cright%29)

### (3) Equation of Motion
![formula](https://latex.codecogs.com/svg.latex?m%5Cfrac%7Bdv%7D%7Bdt%7D%20%3D%20mg%20-%20F_d%20-%20F_%7B%5Ctext%7Blin%7D%7D%20-%20F_%7B%5Ctext%7Bmag%7D%7D)
![formula](https://latex.codecogs.com/svg.latex?F_d%20%3D%20%5Cfrac%7B1%7D%7B2%7D%5Crho%20C_d%20A%20v%7Cv%7C%2C%20%5Cquad%20F_%7B%5Ctext%7Blin%7D%7D%20%3D%20k_%7B%5Ctext%7Blin%7D%7Dv)

### (4) RL Circuit & Implicit Integration
![formula](https://latex.codecogs.com/svg.latex?L%5Cfrac%7Bdi%7D%7Bdt%7D%20%2B%20Ri%20%3D%20%5Ctext%7BEMF%7D)
![formula](https://latex.codecogs.com/svg.latex?i_k%20%3D%20%5Cfrac%7Bi_%7Bk-1%7D%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7BL%7D%5Ctext%7BEMF%7D_%7Bk-1%7D%7D%7B1%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7BL%7DR_%7B%5Ctext%7Btot%7D%7D%7D)

### (5) Lenz Force (Power balance)
![formula](https://latex.codecogs.com/svg.latex?F_%7B%5Ctext%7Bmag%7D%7D%20%3D%20%5Cfrac%7B%5Ctext%7BEMF%7D%5Ccdot%20i%7D%7Bv%7D)

### (6) Wheeler Approximation (Inductance)
![formula](https://latex.codecogs.com/svg.latex?L%28%5Cmu%20H%29%20%5Capprox%20%5Cfrac%7BN%5E2%20D%5E2%7D%7B18D%20%2B%2040%5Cell%7D)

---

## 3. Installation & Execution

```bash
pip install streamlit numpy pandas matplotlib
streamlit run emf_lab_app_v2.py
```

---

## 4. Program Arguments

| **Category** | **Argument** | **Unit** | **Default** | **Description** |
|---------------|--------------|-----------|--------------|-----------------|
| **Simulation** | `--h0` | m | 0.9 | Drop height (tube length) |
|  | `--z0` | m | 0.0 | Initial position |
|  | `--v0` | m/s | 0.0 | Initial velocity |
|  | `--t_max` | s | 1.5 | Total simulation time |
|  | `--dt` | s | 5e-4 | Time step size |
|  | `--dt_sub` | - | 5 | Substeps for RL integration |
| **Magnet** | `--mag_m` | kg | 0.025 | Magnet mass |
|  | `--mag_R` | m | 0.006 | Magnet radius |
|  | `--mag_L` | m | 0.010 | Magnet length |
|  | `--Br` | T | 1.2 | Residual flux density |
| **Environment** | `--rho` | kg/m³ | 1.2 | Air density |
|  | `--cd` | - | 0.6 | Drag coefficient |
|  | `--k_lin` | N·s/m | 0.0 | Linear damping |
|  | `--tube_R` | m | 0.012 | Tube radius |
| **Circuit** | `--mode` | - | `series` | Circuit type (`open`/`series`) |
|  | `--R_load` | Ω | 50 | Load resistance |
|  | `--voltage_unit` | - | `mV` | Voltage display unit |
|  | `--results_dir` | - | `results` | Result folder |
| **Coil** | `--coil` | - | - | Coil definition string (multiple allowed) |

---

## 5. Output Files
- `lab_run_timeseries.csv` — time series (z, v, EMF, current, voltage, Lenz force)
- `lab_run_peaks.csv` — peak summary table
- `lab_run_plot_emf.png` — induced EMF
- `lab_run_plot_v.png` — velocity vs time
- `lab_run_plot_Fmag.png` — Lenz force
- `lab_run_plot_z.png` — position vs time

---

## 6. References
| Property | Symbol | Value | Source |
|-----------|---------|--------|---------|
| Vacuum permeability | μ₀ | 4π×10⁻⁷ H/m | NIST CODATA (2019) |
| Copper resistivity | ρ(Cu) | 1.68×10⁻⁸ Ω·m | NIST Material Data |
| Air density | ρ | 1.2 kg/m³ | 20 °C, sea level |
| Gravity | g | 9.80665 m/s² | International Standard |
| Drag coefficient | C_d | 0.5–0.7 | NASA Glenn Data |
| NdFeB residual flux | B_r | 1.0–1.48 T | Stanford Magnets / Arnold |
| Temperature coefficient | α_T | –0.08 ~ –0.12 %/°C | Stanford Magnets |
| Inductance formula | L | Wheeler approximation | Wheeler (1928) |

---

> © 2025 Aiden Jung 정윤서 — Seoul Science High School
> *EMF Lab Simulator v4.1 — Numerical Validation of Faraday’s and Lenz’s Laws*