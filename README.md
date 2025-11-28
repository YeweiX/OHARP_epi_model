ESBL Transmission Dynamics and Burden Model
Authors: Yewei Xie
Affiliation: Duke-NUS Medical School
Date: November 2025

Overview
This repository contains the R and Stan code used to simulate the transmission dynamics of Extended-Spectrum Beta-Lactamase (ESBL) producing Enterobacteriaceae. The analysis is performed using a compartmental ODE model calibrated to longitudinal surveillance data.
The repository allows reviewers to reproduce:
Bayesian Calibration: Fitting the model to observed incidence data.
Scenario Analysis: Estimating the burden of resistance via counterfactual simulations.
Sensitivity Analysis: Assessing parameter uncertainty using One-Way (OWSA) and Global (PRCC) methods.
Repository Structure
The scripts are numbered to indicate the order of execution.
code
Text
.
├── model/
│   └── esbl_model.stan             # Stan model file (Bayesian ODE system)
├── scripts/
│   ├── 01_calibration.R            # Main calibration script (runs MCMC)
│   ├── 02_burden_analysis.R        # Counterfactual analysis (Baseline vs. No-Resistance)
│   ├── 03_sensitivity_owsa.R       # One-Way Sensitivity Analysis (Parameter sweeps)
│   └── 04_sensitivity_prcc.R       # Global Sensitivity Analysis (LHS-PRCC)
├── outputs/                        # Directory for saved model objects and plots
└── README.md                       # This file
System Requirements
Software
R (v4.0.0 or higher)
RStudio (Recommended)
C++ Compiler (Required for rstan):
Windows: Rtools (version matching your R installation).
macOS: Xcode Command Line Tools.
Linux: g++ or clang++.
Required Packages
Run the following R command to install the necessary dependencies:
code
R
install.packages(c("deSolve", "rstan", "dplyr", "tidyr", "ggplot2", 
                   "patchwork", "scales", "lhs", "sensitivity", 
                   "doParallel", "foreach", "gridExtra"))
Note: For rstan, ensure your C++ toolchain is correctly configured.
Instructions for Reproduction
To reproduce the model results, please execute the scripts in the following sequence. Ensure your working directory is set to the root of this repository.
1. Model Calibration
Script: scripts/01_calibration.R
Action: This script compiles the Stan model, loads the surveillance data, and performs Hamiltonian Monte Carlo (HMC) sampling.
Output:
Saves the fitted model object (posterior samples) to outputs/fit_esbl_results.RData.
Generates diagnostic trace plots and posterior predictive checks.
Note: This step is computationally intensive and may take 30–60 minutes depending on hardware.
2. Burden Estimation (Counterfactuals)
Script: scripts/02_burden_analysis.R
Prerequisite: Requires fit_esbl_results.RData generated in Step 1.
Action: Simulates two parallel scenarios (Baseline vs. No-Resistance) over a 10-year horizon using the posterior parameter sets. It calculates the incremental mortality and morbidity attributable to resistance.
Output: Generates time-series plots of averted deaths and summary statistics for hospital and community settings.
3. One-Way Sensitivity Analysis
Script: scripts/03_sensitivity_owsa.R
Action: Identifies the top influential parameters via a pre-scan (Tornado method) and performs a detailed parameter sweep (+/- 50%) for the most influential parameters.
Output: Generates spider plots showing the impact of individual parameter variations on cumulative outcomes.
4. Global Sensitivity Analysis (LHS-PRCC)
Script: scripts/04_sensitivity_prcc.R
Action: Performs Latin Hypercube Sampling (LHS) to generate 1,000 parameter sets, runs parallel simulations, and calculates Partial Rank Correlation Coefficients (PRCC).
Output: Generates a heatmap visualizing the correlation between biological parameters and key model outcomes (cumulative infections and deaths).
Data Availability
Surveillance Data: Aggregate incidence data used for calibration is hardcoded within 01_calibration.R to ensure immediate reproducibility. No individual patient-level data is included.
Demographics: Parameters for birth, death, and migration rates are derived from publicly available national statistics (sources cited in the manuscript).
Contact
For technical queries regarding code execution, please contact:
yewei.xie@u.duke.nus.edu
