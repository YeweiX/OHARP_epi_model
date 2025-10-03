# OHARP_epi_model
Epidemiological Modeling of ESBL Infection Burden
This project contains R scripts to model the burden of Extended-Spectrum Beta-Lactamase (ESBL) producing bacterial infections, including model generation, posterior checks, scenario analysis, and sensitivity analyses.

Project Structure
The project is organized into several R scripts, each responsible for a specific part of the analysis pipeline. The core of the modeling is performed using Stan for MCMC simulation.

File Descriptions
stancode_opencohort_norisk_4.stan: The Stan code defining the statistical model.

model_mcmc_0906.R: R script to run the MCMC simulation using the Stan model to generate model parameters.

modelposterior.R: R script for conducting posterior predictive checks to evaluate the model fit.

epi_scenario0708_v2.R: R script to calculate and generate the baseline epidemiological burden based on the model outputs.

epi_scenario_overall_0829.R: R script to calculate the incremental burden of ESBL infections under different scenarios.

osa_epi_0724.R: R script to perform a one-way sensitivity analysis (OSA) on the epidemiological model.

fig5.R: R script used to generate plots for the one-way sensitivity analysis results.

prcc_0704_v2.R: R script to perform a Partial Rank Correlation Coefficient (PRCC) analysis and generate corresponding plots.

How to Run the Analysis
Follow these steps in order to reproduce the analysis.

1. MCMC Model Generation
First, run the main model to generate the MCMC samples. This script utilizes the Stan code to fit the model to the data.

Rscript model_mcmc_0906.R

2. Model Posterior Check
After running the model, perform a posterior check to ensure the model has converged and is a good fit for the data.

Rscript modelposterior.R

3. Burden Calculation
Next, calculate the baseline and incremental burden of ESBL infections.

Baseline Burden:

Rscript epi_scenario0708_v2.R

Incremental Burden:

Rscript epi_scenario_overall_0829.R

4. Sensitivity Analyses
Finally, run the sensitivity analyses to understand the impact of different parameters on the model outputs.

One-Way Sensitivity Analysis (OSA):
First, run the analysis script:

Rscript osa_epi_0724.R

Then, generate the corresponding plot:

Rscript fig5.R

PRCC Analysis:
Run the PRCC analysis and generate the plot:

Rscript prcc_0704_v2.R

Dependencies
This analysis is conducted in R. You will need to have R installed, along with several packages. Key dependencies include:

rstan

ggplot2

dplyr

(Please add any other packages used in the scripts here)

You can install the necessary packages in R using install.packages("package_name").
