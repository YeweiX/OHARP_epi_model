## Overview

This repository contains source code to reproduce a mathematical modelling analysis of Extended-Spectrum Beta-Lactamase (ESBL)-producing Enterobacteriaceae transmission dynamics and burden estimation.

The analysis uses a dynamic compartmental ODE model to calibrate transmission parameters to longitudinal surveillance data using Bayesian inference (HMC via Stan), estimate burden via counterfactual simulations (Baseline versus No-Resistance), and assess uncertainty using both one-way sensitivity analysis (OWSA) and global sensitivity analysis (LHS-PRCC).

## Repository Structure

Scripts are numbered in the intended execution order.

```text
.
├── model/
│   └── esbl_model.stan             # Stan model file (ODE system & likelihood)
├── scripts/
│   ├── 01_calibration.R            # Model calibration (MCMC sampling)
│   ├── 02_burden_analysis.R        # Counterfactual burden analysis
│   ├── 03_sensitivity_owsa.R       # One-way sensitivity analysis (spider plots)
│   └── 04_sensitivity_prcc.R       # Global sensitivity analysis (LHS-PRCC heatmap)
├── outputs/                        # Saved model objects (.RData) and derived outputs
└── README.md                       # This file
````

## System Requirements

### Software

* **R** (>= 4.0.0 recommended)
* **RStudio** (optional, recommended)
* **C++ toolchain** (required to compile Stan models):

  * **Windows:** Rtools (matching your R version)
  * **macOS:** Xcode Command Line Tools


### R Dependencies

Install required packages:

```r
install.packages(c(
  "deSolve", "rstan", "dplyr", "tidyr", "ggplot2",
  "patchwork", "scales", "lhs", "sensitivity",
  "doParallel", "foreach", "gridExtra"
))
```

> Note: Ensure `rstan` is correctly configured with your C++ compiler before running `scripts/01_calibration.R`.

## Reproduction Workflow

Run scripts in numeric order. Set your working directory to the repository root.

### 1) Model Calibration

**File:** `scripts/01_calibration.R`

**What it does:** Compiles the Stan model, loads aggregate surveillance data, and runs Hamiltonian Monte Carlo (HMC) to estimate posterior parameter distributions.

**Outputs:**

* Saves fitted model object to `outputs/fit_esbl_results.RData`
* Generates diagnostic plots (trace plots, posterior predictive checks)

### 2) Burden Estimation (Counterfactuals)

**File:** `scripts/02_burden_analysis.R`
**Prerequisite:** Requires `outputs/fit_esbl_results.RData` from Step 1.

**What it does:** Uses posterior samples to simulate two scenarios over a 10-year horizon:

* **Baseline:** calibrated ESBL transmission dynamics
* **Counterfactual:** theoretical **No-Resistance** scenario

**Outputs:** Incremental burden estimates (e.g., net attributable mortality/morbidity) and summary visualizations.

### 3) One-Way Sensitivity Analysis (OWSA)

**File:** `scripts/03_sensitivity_owsa.R`

**What it does:** Performs a tornado-style pre-scan to identify influential parameters, then sweeps key parameters (e.g., ±50%) to visualize outcome elasticity.

**Outputs:** Spider plots for key outcomes (e.g., hospital deaths, community cases).

### 4) Global Sensitivity Analysis (LHS-PRCC)

**File:** `scripts/04_sensitivity_prcc.R`

**What it does:** Generates parameter sets using Latin Hypercube Sampling (e.g., 1,000 draws) and computes Partial Rank Correlation Coefficients (PRCC), using parallel simulation where configured.

**Outputs:** Heatmap ranking correlations between input parameters and cumulative 10-year outcomes.

## Data Availability

* **Surveillance data:** Aggregate hospital incidence data used for calibration is embedded directly in `scripts/01_calibration.R` for reproducibility. No individual patient-level data is included.
* **Demographics:** Birth, death, and migration parameters come from publicly available national statistics (see manuscript for details).

## Contact
For technical questions about the code or reproducibility, please contact the corresponding author listed in the manuscript.

