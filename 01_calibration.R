# ==============================================================================
# Analysis Script: ESBL Transmission Model (OHARP)
# Description: Bayesian calibration of a compartmental transmission model using Stan.
# Date: November 2025
# ==============================================================================

# --- 1. Setup and Environment -------------------------------------------------

# Clear workspace
rm(list = ls())

# Set Working Directory
# NOTE: Set this to the folder containing your .stan file and data
working_dir <- "." 
setwd(working_dir)

# Load required packages
required_packages <- c("deSolve", "ggplot2", "tidyr", "gridExtra", 
                       "StanHeaders", "rstan", "bayesplot", "dplyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# Configure Stan options for parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# --- 2. Demographic and Initial Conditions (2011 Baseline) --------------------

# Population sizes
N_total <- 5.18e6          # Total resident population (Singapore, 2011)
N_hosp  <- 37588 * 0.857   # Estimated active hospital population
N_comm  <- N_total - N_hosp # Community population

# --- Compartment Initialization ---
# Notation: 
# S = Susceptible, RC = Resistant Colonized, RI = Resistant Infected, SI = Sensitive Infected
# Suffix: _h = Hospital, _c = Community

# Hospital Compartments (derived from prevalence data)
RC_h <- round(0.124 * N_hosp)
RI_h <- round(0.119 * 0.314 * 0.377 * N_hosp)
SI_h <- round(0.119 * 0.314 * 0.560 * N_hosp)
S_h  <- N_hosp - RC_h - RI_h - SI_h

# Community Compartments (derived from prevalence data)
RC_c <- round(0.257 * N_comm)
RI_c <- round(0.001 * N_comm)
SI_c <- round(0.010 * N_comm)
S_c  <- N_comm - RC_c - RI_c - SI_c

# Cumulative counters (Initial states = 0)
# Dr = Death due to resistance, Ds = Death due to sensitive, C_ = Cumulative incidence
Dr_h <- 0; Ds_h <- 0
Dr_c <- 0; Ds_c <- 0
Cum_RI_h <- 0; Cum_SI_h <- 0
Cum_RI_c <- 0; Cum_SI_c <- 0

# Vector of Initial Conditions (y0)
y0 <- c(
  Sh    = S_h,
  RCh   = RC_h,
  RIh   = RI_h,
  SIh   = SI_h,
  Sc    = S_c,
  RCc   = RC_c,
  RIc   = RI_c,
  SIc   = SI_c,
  Drh   = Dr_h,
  Dsh   = Ds_h,
  Drc   = Dr_c,
  Dsc   = Ds_c,
  C_RIh = Cum_RI_h,
  C_SIh = Cum_SI_h,
  C_SIc = Cum_SI_c,
  C_RIc = Cum_RI_c
)

# Consistency Check
if(abs(sum(y0[1:12]) - N_total) > 1) {
  warning("Initial conditions do not sum to total population.")
} else {
  message("Initial conditions verified against total population.")
}


# --- 3. Observed Data (2012–2021) ---------------------------------------------

# Incidence per 10,000 inpatient days
incidence_raw <- c(38, 35, 34.5, 33.5, 37, 27.5, 28.5, 28.2, 27, 24.5)

obs_data <- data.frame(
  year = 2012:2021,
  RIh_obs = round(incidence_raw)
)

RIh_obs <- obs_data$RIh_obs
n_years <- nrow(obs_data)
n_days  <- 365 * n_years
t0      <- 0
ts      <- seq(1, n_days, by = 1)


# --- 4. Model Parameters ------------------------------------------------------

# 4.1. Hospital Parameters
mu_sih       <- 0.00311667  # Mortality rate (Sensitive Infection, Hospital)
mu_rih       <- 0.00404500  # Mortality rate (Resistant Infection, Hospital)
beta_sh_sih  <- 0.00348     # Transmission rate (S -> SI, Hospital)

# Hospital Recovery/Transition Rates (gamma)
gamma_sih_sh  <- 0.142857   # Recovery: SI -> S
gamma_rch_sh  <- 0.00612    # Clearance: RC -> S
gamma_rih_rch <- 0.1111     # Recovery: RI -> RC

# 4.2. Community Parameters
beta_sc_sic   <- 0.00000792 # Transmission rate (S -> SI, Community)

# Community Recovery/Transition Rates
gamma_sic_sc  <- 0.25       # Recovery: SI -> S
gamma_rcc_sc  <- 0.009      # Clearance: RC -> S
gamma_ric_rcc <- 0.17       # Recovery: RI -> RC

# Community Mortality
mu_sic <- 0.003226
mu_ric <- 0.003629

# 4.3. Movement & Progression
alpha_adm     <- 0.000316     # Hospital Admission Rate
alpha_dis     <- 0.235        # Hospital Discharge Rate
delta_rch_sih <- 0.00084658   # Progression: RC -> SI (Hospital)
delta_rcc_sic <- 0.00084658   # Progression: RC -> SI (Community)

# 4.4. Demographics (Daily Rates, SingStat 2011-2021)
mu_b <- 0.00002518  # Crude Birth Rate
mu_m <- 0.000002    # Net Migration Rate
mu_d <- 0.00001353  # Crude Death Rate


# --- 5. Bayesian Inference (Stan) ---------------------------------------------

# Prepare Data List
stan_data <- list(
  n_years = n_years,
  n_days  = n_days,
  y0      = y0,
  t0      = t0,
  ts      = ts,
  RIh_obs = RIh_obs,
  
  # Fixed Parameter Priors
  mu_sih_fixed        = mu_sih,
  mu_rih_fixed        = mu_rih,
  gamma_sih_sh_fixed  = gamma_sih_sh,
  gamma_rch_sh_fixed  = gamma_rch_sh,
  gamma_rih_rch_fixed = gamma_rih_rch,
  beta_sh_sih_fixed   = beta_sh_sih,
  beta_sc_sic_fixed   = beta_sc_sic,
  mu_sic_fixed        = mu_sic,
  mu_ric_fixed        = mu_ric,
  gamma_sic_sc_fixed  = gamma_sic_sc,
  gamma_rcc_sc_fixed  = gamma_rcc_sc,
  gamma_ric_rcc_fixed = gamma_ric_rcc,
  delta_rcc_sic_fixed = delta_rcc_sic,
  delta_rch_sih_fixed = delta_rch_sih,
  alpha_adm_fixed     = alpha_adm,
  alpha_dis_fixed     = alpha_dis,
  mu_b_fixed          = mu_b,
  mu_m_fixed          = mu_m,
  mu_d_fixed          = mu_d
)

# Define Stan Model Path (Relative)
stan_file_path <- "stancode_opencohort_norisk_4.stan"

# Check if Stan file exists
if (!file.exists(stan_file_path)) {
  stop("Stan file not found in working directory: ", stan_file_path)
}

# Run MCMC Sampling
message("Starting Stan model calibration...")
fit_esbl <- stan(
  file    = stan_file_path,
  data    = stan_data,
  iter    = 4000,
  chains  = 4,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

# Save the fitted model
save(fit_esbl, file = "fit_esbl_results.RData")


# --- 6. Post-Processing and Diagnostics ---------------------------------------

# Summary of Parameters
params_of_interest <- c("beta_hh", "beta_eh", "beta_hc", "beta_sc_ac", 
                        "beta_sc_ec", "beta_sh_sih", "delta_rcc_ric")

print(fit_esbl, pars = params_of_interest, digits = 4)

# Diagnostic Plots
trace_plot <- traceplot(fit_esbl, pars = params_of_interest)
print(trace_plot)

posterior_density <- stan_dens(fit_esbl, pars = params_of_interest, separate_chains = TRUE)
print(posterior_density)


# --- 7. Visualization: Posterior Predictive Check -----------------------------

# Extract Posterior Draws
# Note: Ensure 'pred_yearly_total_RIh_incidence' matches the variable name in your .stan file
incidence_post <- as.matrix(fit_esbl, pars = "pred_yearly_total_RIh_incidence")

# Set Color Scheme
color_scheme_set("brightblue")

# Generate Posterior Predictive Plot
ppc_plot <- ppc_ribbon(
  y = obs_data$RIh_obs,            # Observed data
  yrep = incidence_post,           # Simulated data (posterior predictive)
  x = 1:nrow(obs_data),            # X-axis (Year Index)
  prob = 0.5,                      # 50% Interval (inner)
  prob_outer = 0.95                # 95% Interval (outer)
) +
  geom_point(
    aes(x = 1:nrow(obs_data), y = obs_data$RIh_obs), 
    color = "navy", 
    size = 3
  ) +
  scale_x_continuous(
    breaks = 1:nrow(obs_data), 
    labels = obs_data$year
  ) +
  labs(
    x = "Year",
    y = "Incidence per 10,000 inpatient days",
    title = "Model Calibration: Observed vs. Predicted Incidence",
    subtitle = "Blue points: Observed Data | Ribbon: 50% and 95% Posterior Credible Intervals"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  )

# Display and Save Plot
print(ppc_plot)
ggsave("Figure_Calibration_PPC.png", plot = ppc_plot, width = 10, height = 6, dpi = 300)