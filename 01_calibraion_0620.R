# Create a new environment
esbl_env <- new.env()

# Attach the environment to work within it
attach(esbl_env)
setwd("C:/Users/e0977434/Desktop/OneDrive/PhD thesis/Models/OHARP model/Codes/code0523")
# Install packages if not already installed
if (!requireNamespace("deSolve", quietly = TRUE)) install.packages("deSolve")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("StanHeaders", quietly = TRUE)) install.packages("StanHeaders")
if (!requireNamespace("rstan", quietly = TRUE)) install.packages("rstan")

# Configure rstan
library(deSolve)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(StanHeaders)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# --- 1. Define Parameters ---

# Population & Model Initial Conditions (2011 for 2012 Start)
N <- 5.18e6         # Total resident population (Singapore, 2011)
Nh <- 37588*0.857         # Hospital population (fixed)
Nc <- N - Nh        # Community population

# Hospital Initial Conditions (Consolidated)
RCh <- round(0.124 * Nh)
RIh <- round(0.119 * 0.314 * 0.377 * Nh)
SIh <- round(0.119 * 0.314 * 0.56 * Nh)
Sh <- Nh - RCh - RIh - SIh
Drh <- 0 # Initial cumulative deaths are 0
Dsh <- 0 # Initial cumulative deaths are 0

# Community Initial Conditions (Consolidated)
RCc <- round(0.257 * Nc)
RIc <- round(0.001 * Nc)
SIc <- round(0.01 * Nc)
Sc <- Nc - RCc - RIc - SIc
Drc <- 0 # Initial cumulative deaths are 0
Dsc <- 0 # Initial cumulative deaths are 0

# Cumulative resistant infections (hospital) - only one needed for calibration
C_RIh <- 0
C_SIh <- 0 # Cumulative resistant colonisations (hospital) - only one needed for calibration

# Initial conditions vector (13 compartments - matching new Stan model)
y0 <- c(
  Sh = Sh,
  RCh = RCh,
  RIh = RIh,
  SIh = SIh,
  Sc = Sc,
  RCc = RCc,
  RIc = RIc,
  SIc = SIc,
  Drh = Drh,
  Dsh = Dsh,
  Drc = Drc,
  Dsc = Dsc,
  C_RIh = C_RIh,
  C_SIh = C_SIh,
  C_SIc = 0,  # Cumulative resistant colonisations (community) - only one needed for calibration
  C_RIc = 0   # Cumulative resistant infections (community) - only one needed for calibration
)

# -----------------------------
# Print revised initial conditions:
print(y0)

total <- Sh + RCh + RIh + SIh + Drh + Dsh+ Sc + RCc + RIc + SIc + Drc + Dsc
if(total != N) {
  warning("The sum of compartments (", total, ") does not equal the total population (", N, ").")
} else {
  print("Initial conditions are consistent with the total population.")
}


# -----------------------------
# Observed Data (2012–2021)
# -----------------------------
incidence_per_10k <- c(38, 35, 34.5, 33.5, 37, 27.5, 28.5, 28.2, 27, 24.5)
yearly_data <- data.frame(
  year = 2012:2021,
  RIh_obs = round(incidence_per_10k)
)
RIh_obs <- round(incidence_per_10k)
n_years <- nrow(yearly_data)  # 10 years
n_days <- 365 * n_years       # 3650 days (2012–2021)
t0 <- 0
ts <- seq(1, n_days, by = 1)



# -----------------------------
# Parameters
# -----------------------------
# Hospital Parameters
mu_sih <- 0.00311667
mu_rih <- 0.00404500
delta_rch_rih <- 0.03395087

# Hospital Recovery Rates
gamma_sih_sh <- 0.142857
gamma_rch_sh <- 0.00612
gamma_rih_rch <- 0.1111

# Community Parameters
beta_sc_sic <- 0.00000792  


# Community Recovery Rates
gamma_sic_sc <- 0.25
gamma_rcc_sc <- 0.009
gamma_ric_rcc <- 0.17

# Community Mortality
mu_sic <- 0.003226
mu_ric <- 0.003629

# Transition Probabilities
alpha_adm <- 0.000316
alpha_dis <- 0.235

# Demographic Parameters
mu_b <- 0.00002518    # Average Crude Birth Rate (Daily) 2011-2021, derived from SingStat
mu_m <- 0.000002      # Average Net Migration Rate (Daily) 2011-2021, derived from SingStat
mu_d <- 0.00001353      # Average Crude Death Rate (Daily) 2011-2021, derived from SingStat

delta_rch_sih<-0.00084658
delta_rcc_sic<-0.00084658



# -----------------------------
# Prepare Data for Stan
# -----------------------------
data_sir <- list(
  n_years = n_years,
  n_days = n_days,
  y0 = y0,
  t0 = t0,
  ts = ts,
  RIh_obs = RIh_obs,

  # Fixed parameters
  mu_sih_fixed = mu_sih,
  mu_rih_fixed = mu_rih,
  gamma_sih_sh_fixed = gamma_sih_sh,
  gamma_rch_sh_fixed = gamma_rch_sh,
  gamma_rih_rch_fixed = gamma_rih_rch,
  beta_sc_sic_fixed = beta_sc_sic,
  mu_sic_fixed = mu_sic,
  mu_ric_fixed = mu_ric,
  gamma_sic_sc_fixed = gamma_sic_sc,
  gamma_rcc_sc_fixed = gamma_rcc_sc,
  gamma_ric_rcc_fixed = gamma_ric_rcc,
  delta_rcc_sic_fixed = delta_rcc_sic,
  delta_rch_sih_fixed = delta_rch_sih,
  delta_rch_rih_fixed = delta_rch_rih,
  alpha_adm_fixed = alpha_adm,
  alpha_dis_fixed = alpha_dis,
  mu_b_fixed = mu_b,
  mu_m_fixed = mu_m,
  mu_d_fixed = mu_d
)

# Print Stan data
print(data_sir)

# -----------------------------
# Run Stan Model
# -----------------------------
options(mc.cores = parallel::detectCores())

fit6 <- stan(
  file = "C:/Users/e0977434/Desktop/OneDrive/PhD thesis/Models/OHARP model/Codes/stancode_opencohort_norisk_0923.stan",  # for windows
  data = data_sir,
  iter = 4000,
  chains = 4,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

save(fit6, file = "fit60922.RData")  # Save the fitted model object for later use

 

fit5 <- stan(
  file = "/Users/yeweixie/Library/CloudStorage/OneDrive-Personal(2)/PhD thesis/Models/OHARP model/Codes/stancode_opencohort_norisk_3.stan",  # for mac
  data = data_sir,
  iter = 1000,
  chains = 4,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

# -----------------------------
# Post-Processing

load("/Users/yeweixie/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/fit60922.RData")  # Load the fitted model object if needed
# -----------------------------
# Check convergence and diagnostics
pars1 = c("beta_hh", "beta_eh", "beta_hc", "beta_sc_ac", "beta_sc_ec", "beta_sh_sih", "delta_rcc_ric")



##Provides a quick summary of the parameters and diagnostics
traceplot(fit6, pars=pars1)
stan_dens(fit6, pars=pars1, separate_chains=TRUE)
print(fit6, digits = 9, pars=pars1)  #
save(fit6, file = "fit6.RData")  # Save the fitted model object for later use

draws <- as.data.frame(fit6, pars = pars1)  # each column is a parameter
minmax <- t(vapply(draws, function(x) c(min = min(x), max = max(x)), c(min=0, max=0)))
print(round(minmax, 9))
###

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(rstanarm)

# --- Data Preparation ---
# Extract posterior draws for the simulated incidence
# Assuming 'fit5' is your fitted model object and 'sim_RIh' is the relevant parameter
incidence_post <- as.matrix(fit6, pars = "pred_yearly_total_RIh_incidence")
yrep<-incidence_post

R <- as.matrix(fit6, pars = "average_community_resistant_proportion")
hair<- as.matrix(fit6, pars = "average_total_hai_sih_ratio")


# Get the observed incidence data
# Assuming 'yearly_data' is your data frame and 'RIh_obs' is the column with observed values
obs_incidence <- yearly_data$RIh_obs

# Create a data frame for observed points to facilitate legend mapping
obs_df <- data.frame(
  Year = 1:length(obs_incidence),
  Incidence = obs_incidence,
  DataType = "Observed (y)" # Add a column to map to color aesthetic
)

# Assuming 'obs_incidence' and 'incidence_post' are loaded
color_scheme_set("brightblue")
# --- Plotting ---
ppc_ribbon(
  y = obs_incidence,           # Provide the actual observed data
  yrep = incidence_post,       # Provide the matrix of posterior predictive simulations
  x = 1:length(obs_incidence), # Use observation index (Year) for x-axis
  # y_draw = "line",           # Optional: tell ppc_ribbon to only draw the line for y
  # size = 0.5,                # Optional: Adjust size for elements drawn by ppc_ribbon (e.g., ribbon lines) if needed
  prob = 0.5,                  # Probability for the thicker inner intervals
  prob_outer = 0.95             # Probability for the thinner outer intervals
) +
  # Add observed data points manually with custom appearance
  geom_point(
    aes(x = 1:length(obs_incidence), y = obs_incidence), # Map x and y to observed data
    color = "navy",  # Set desired color for points
    size = 3        # Set desired size for points (adjust as needed)
  ) +
  # Add labels and theme elements
  labs(
    x = "Year",
    y = "Incidence per 10,000 inpatient days",
    title = "Posterior predictive intervals vs Observed Incidence",
    subtitle = "Showing 50% & 95% intervals (yrep) and Observed data (y, blue points)" # Updated subtitle
  ) +
  theme_minimal() # Or your preferred theme
# Save the plot
ggsave("ppc_ribbon_plot.png", width = 10, height = 6)
# -----------------------------

