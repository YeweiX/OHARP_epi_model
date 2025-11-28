# ==============================================================================
# Analysis: Global Sensitivity Analysis (LHS-PRCC)
# Description: Evaluates the global influence of uncertainty in all parameters
#              simultaneously on key model outcomes using Partial Rank 
#              Correlation Coefficients (PRCC).
# Methodology: Latin Hypercube Sampling (LHS) -> Parallel Simulation -> PRCC
# ==============================================================================

# --- 1. SETUP AND INITIALIZATION ----------------------------------------------
rm(list = ls())

# Load Libraries
library(deSolve)      # ODE Solver
library(lhs)          # Latin Hypercube Sampling
library(sensitivity)  # PRCC Calculation
library(dplyr)        # Data Manipulation
library(ggplot2)      # Visualization
library(doParallel)   # Parallel Processing
library(foreach)      # Parallel Loops
library(gridExtra)    # Plot Arrangement

# Parallel Setup
num_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Parallel processing initialized on", num_cores, "cores.\n")
set.seed(123)

# --- 2. MODEL DEFINITION ------------------------------------------------------

# 2.1 Initial Conditions (Baseline 2011)
N_total <- 5.18e6
N_hosp  <- 37588 * 0.857
N_comm  <- N_total - N_hosp

# Hospital Compartments
RCh <- round(0.124 * N_hosp)
RIh <- round(0.119 * 0.314 * 0.377 * N_hosp)
SIh <- round(0.119 * 0.314 * 0.56 * N_hosp)
Sh  <- N_hosp - RCh - RIh - SIh

# Community Compartments
RCc <- round(0.257 * N_comm)
RIc <- round(0.001 * N_comm)
SIc <- round(0.01 * N_comm)
Sc  <- N_comm - RCc - RIc - SIc

y_init <- c(
  Sh = Sh, RCh = RCh, RIh = RIh, SIh = SIh,
  Sc = Sc, RCc = RCc, RIc = RIc, SIc = SIc,
  Drh = 0, Dsh = 0, Drc = 0, Dsc = 0,
  C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0
)

# 2.2 ODE System
constrain_flow <- function(rate, pool) { min(rate, pool) }

infection_model_ode <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Populations
    Nh_curr <- max(1, Sh + RCh + RIh + SIh)
    Nc_curr <- max(1, Sc + RCc + RIc + SIc)
    N_living <- Nh_curr + Nc_curr
    
    # --- Hospital Flows ---
    f_Sh_env <- constrain_flow(beta_eh * Sh, Sh)
    f_Sh_rc  <- constrain_flow(beta_hh * Sh * RCh / Nh_curr, Sh)
    f_Sh_si  <- constrain_flow(beta_sh_sih * Sh * SIh / Nh_curr, Sh)
    f_Sh_out <- constrain_flow(alpha_dis * Sh, Sh)
    f_Sh_d   <- constrain_flow(mu_d * Sh, Sh)
    
    f_RCh_ri  <- constrain_flow(delta_rch_rih * RCh, RCh)
    f_RCh_s   <- constrain_flow(gamma_rch_sh * RCh, RCh)
    f_RCh_out <- constrain_flow(alpha_dis * RCh, RCh)
    f_RCh_si  <- constrain_flow(delta_rch_sih * RCh, RCh)
    f_RCh_d   <- constrain_flow(mu_d * RCh, RCh)
    
    f_RIh_rc <- constrain_flow(gamma_rih_rch * RIh, RIh)
    f_RIh_d  <- constrain_flow(mu_rih * RIh, RIh)
    
    f_SIh_s  <- constrain_flow(gamma_sih_sh * SIh, SIh)
    f_SIh_d  <- constrain_flow(mu_sih * SIh, SIh)
    
    # --- Community Flows ---
    f_Sc_env   <- constrain_flow(beta_sc_ec * Sc, Sc)
    f_Sc_ani   <- constrain_flow(beta_sc_ac * Sc * 0.33, Sc)
    f_Sc_rc    <- constrain_flow(beta_hc * Sc * RCc / Nc_curr, Sc)
    f_Sc_si    <- constrain_flow(beta_sc_sic * Sc * SIc / Nc_curr, Sc)
    f_Sc_in    <- constrain_flow(alpha_adm * Sc, Sc)
    f_Sc_d     <- constrain_flow(mu_d * Sc, Sc)
    
    f_RCc_ri   <- constrain_flow(delta_rcc_ric * RCc, RCc)
    f_RCc_si   <- constrain_flow(delta_rcc_sic * RCc, RCc)
    f_RCc_s    <- constrain_flow(gamma_rcc_sc * RCc, RCc)
    f_RCc_in   <- constrain_flow(alpha_adm * RCc, RCc)
    f_RCc_d    <- constrain_flow(mu_d * RCc, RCc)
    
    f_RIc_rc <- constrain_flow(gamma_ric_rcc * RIc, RIc)
    f_RIc_d  <- constrain_flow(mu_ric * RIc, RIc)
    f_RIc_in <- constrain_flow(alpha_adm * RIc, RIc)
    
    f_SIc_s  <- constrain_flow(gamma_sic_sc * SIc, SIc)
    f_SIc_d  <- constrain_flow(mu_sic * SIc, SIc)
    f_SIc_in <- constrain_flow(alpha_adm * SIc, SIc)
    
    # --- Derivatives ---
    dSh  <- -f_Sh_env - f_Sh_rc - f_Sh_si - f_Sh_out - f_Sh_d + f_SIh_s + f_RCh_s + f_Sc_in
    dRCh <- -f_RCh_ri - f_RCh_s - f_RCh_out - f_RCh_si - f_RCh_d + f_Sh_env + f_Sh_rc + f_RIh_rc + f_RCc_in
    dRIh <- -f_RIh_rc - f_RIh_d + f_RCh_ri + f_RIc_in
    dSIh <- -f_SIh_s - f_SIh_d + f_Sh_si + f_RCh_si + f_SIc_in
    
    dSc  <- -f_Sc_env - f_Sc_ani - f_Sc_rc - f_Sc_si - f_Sc_in - f_Sc_d + f_SIc_s + f_RCc_s + f_Sh_out + mu_b * N_living + mu_m * N_living
    dRCc <- -f_RCc_ri - f_RCc_si - f_RCc_s - f_RCc_in - f_RCc_d + f_Sc_env + f_Sc_ani + f_Sc_rc + f_RIc_rc + f_RCh_out
    dRIc <- -f_RIc_rc - f_RIc_d - f_RIc_in + f_RCc_ri
    dSIc <- -f_SIc_s - f_SIc_d - f_SIc_in + f_Sc_si + f_RCc_si
    
    # Cumulative Trackers
    dDrh <- f_RIh_d; dDsh <- f_SIh_d
    dDrc <- f_RIc_d; dDsc <- f_SIc_d
    dC_RIh <- f_RCh_ri + f_RIc_in
    dC_SIh <- f_Sh_si + f_RCh_si + f_SIc_in
    dC_RIc <- f_RCc_ri
    dC_SIc <- f_Sc_si + f_RCc_si
    
    list(c(dSh, dRCh, dRIh, dSIh, dSc, dRCc, dRIc, dSIc,
           dDrh, dDsh, dDrc, dDsc, dC_RIh, dC_SIh, dC_RIc, dC_SIc))
  })
}


# --- 3. PARAMETER SAMPLING (LHS) ----------------------------------------------

# 3.1 Define Parameter Ranges (Uniform Distributions)
param_dists <- list(
  beta_hh       = list(min = 0.002100771, max = 0.013899509),
  beta_eh       = list(min = 0.000000084, max = 0.033214594),
  beta_hc       = list(min = 0.001600901, max = 0.010231657),
  beta_sc_ac    = list(min = 0.000000096, max = 0.003083009),
  beta_sc_ec    = list(min = 0.000000037, max = 0.001111809),
  beta_sh_sih   = list(min = 0.155454327, max = 0.996909995),
  delta_rcc_ric = list(min = 0.000031951, max = 0.000149248),
  delta_rch_rih = list(min = 0.00010513,  max = 0.06779661),
  mu_sih        = list(min = 0.00233333,  max = 0.00390000),
  mu_rih        = list(min = 0.00254333,  max = 0.00554667),
  gamma_sih_sh  = list(min = 0.09090909,  max = 0.20000000),
  gamma_rch_sh  = list(min = 0.00147397,  max = 0.01076667),
  gamma_rih_rch = list(min = 0.06666667,  max = 0.16666667),
  delta_rch_sih = list(min = 0.000592603, max = 0.001100548),
  beta_sc_sic   = list(min = 0.00000707,  max = 0.00000877),
  gamma_sic_sc  = list(min = 0.16666667,  max = 0.33333333),
  gamma_rcc_sc  = list(min = 0.00460000,  max = 0.01800000),
  gamma_ric_rcc = list(min = 0.12500000,  max = 0.25000000),
  mu_sic        = list(min = 0.003226*0.7, max = 0.003226*1.3),
  mu_ric        = list(min = 0.003629*0.7, max = 0.003629*1.3),
  delta_rcc_sic = list(min = 0.000592603, max = 0.001100548),
  alpha_adm     = list(min = 0.000316*0.7, max = 0.000316*1.3),
  alpha_dis     = list(min = 0.235*0.7,    max = 0.235*1.3),
  mu_b          = list(min = 2.518e-5*0.7, max = 2.518e-5*1.3),
  mu_m          = list(min = 2e-6*0.7,     max = 2e-6*1.3),
  mu_d          = list(min = 1.353e-5*0.7, max = 1.353e-5*1.3)
)

# 3.2 Generate LHS Matrix
n_samples <- 1000  # Number of iterations
n_params  <- length(param_dists)
param_names <- names(param_dists)

message("Generating ", n_samples, " parameter sets via LHS...")
lhs_matrix <- randomLHS(n = n_samples, k = n_params)
param_samples <- data.frame(matrix(NA, nrow = n_samples, ncol = n_params))
colnames(param_samples) <- param_names

# Transform Probabilities to Values
for (i in 1:n_params) {
  p_name <- param_names[i]
  p_min  <- param_dists[[p_name]]$min
  p_max  <- param_dists[[p_name]]$max
  param_samples[, i] <- qunif(lhs_matrix[, i], min = p_min, max = p_max)
}

# --- 4. PARALLEL SIMULATION ---------------------------------------------------

# Time Settings
years_total <- 20
years_econ  <- 10
time_steps  <- seq(0, years_total * 365, by = 365)
start_day   <- (years_total - years_econ) * 365

message("Running simulations...")

results_list <- foreach(i = 1:n_samples, .packages = c("deSolve", "dplyr")) %dopar% {
  
  curr_params <- as.list(param_samples[i, ])
  
  # Run ODE
  out <- try(ode(y = y_init, times = time_steps, func = infection_model_ode, 
                 parms = curr_params, method = "lsoda"), silent = TRUE)
  
  if (inherits(out, "try-error")) return(NULL)
  
  df <- as.data.frame(out)
  
  # Calculate Cumulative outcomes over last 10 years
  idx_start <- which(df$time >= start_day)[1]
  if (is.na(idx_start)) return(NULL)
  
  row_start <- df[idx_start, ]
  row_end   <- df[nrow(df), ]
  
  # Return row of outcomes
  data.frame(
    sample_id = i,
    C_RIh = row_end$C_RIh - row_start$C_RIh,
    C_RIc = row_end$C_RIc - row_start$C_RIc,
    Drh   = row_end$Drh   - row_start$Drh,
    Drc   = row_end$Drc   - row_start$Drc
  )
}

stopCluster(cl)
all_results <- bind_rows(results_list)
message("Completed: ", nrow(all_results), "/", n_samples, " runs.")


# --- 5. PRCC CALCULATION ------------------------------------------------------

# Filter Parameters to match successful runs
X <- param_samples[all_results$sample_id, ]
outcomes <- c("C_RIh", "C_RIc", "Drh", "Drc")

# Store PRCC results in a combined dataframe
prcc_list <- list()

for (outcome in outcomes) {
  y <- all_results[[outcome]]
  
  # Compute PRCC (Partial Correlation Coefficients)
  # rank = TRUE makes it PRCC (non-parametric)
  pcc_res <- pcc(X, y, rank = TRUE, nboot = 1000)
  
  temp_df <- pcc_res$PRCC
  temp_df$Parameter <- rownames(temp_df)
  temp_df$Outcome   <- outcome
  prcc_list[[outcome]] <- temp_df
}

all_prcc_df <- bind_rows(prcc_list)

# Save results
write.csv(all_prcc_df, "prcc_results_full.csv", row.names = FALSE)


# --- 6. VISUALIZATION ---------------------------------------------------------

# 6.1 Parameter Label Dictionary (For Pretty Plots)
param_full_names <- c(
  "beta_hh"       = "Hosp. Acquisition (Human)",
  "beta_sh_sih"   = "Hosp. Acquisition (Susceptible)",
  "beta_eh"       = "Hosp. Acquisition (Env)",
  "delta_rch_rih" = "Hosp. Progression (RC->RI)",
  "delta_rch_sih" = "Hosp. Progression (RC->SI)",
  "gamma_sih_sh"  = "Hosp. Recovery (SI->S)",
  "gamma_rih_rch" = "Hosp. Recovery (RI->RC)",
  "gamma_rch_sh"  = "Hosp. Decolonisation",
  "mu_sih"        = "Hosp. Mortality (SI)",
  "mu_rih"        = "Hosp. Mortality (RI)",
  "beta_hc"       = "Comm. Acquisition (Human)",
  "beta_sc_sic"   = "Comm. Acquisition (Susceptible)",
  "beta_sc_ac"    = "Comm. Acquisition (Animal)",
  "beta_sc_ec"    = "Comm. Acquisition (Env)",
  "delta_rcc_ric" = "Comm. Progression (RC->RI)",
  "delta_rcc_sic" = "Comm. Progression (RC->SI)",
  "gamma_sic_sc"  = "Comm. Recovery (SI->S)",
  "gamma_ric_rcc" = "Comm. Recovery (RI->RC)",
  "gamma_rcc_sc"  = "Comm. Decolonisation",
  "mu_sic"        = "Comm. Mortality (SI)",
  "mu_ric"        = "Comm. Mortality (RI)",
  "alpha_adm"     = "Admission Rate",
  "alpha_dis"     = "Discharge Rate",
  "mu_b"          = "Birth Rate",
  "mu_m"          = "Migration Rate",
  "mu_d"          = "Death Rate"
)

# 6.2 Heatmap Preparation
plot_data <- all_prcc_df %>%
  mutate(Parameter_Label = recode(Parameter, !!!param_full_names)) %>%
  filter(original != 0) # Remove empty rows if any

# 6.3 Generate Heatmap
heatmap_plot <- ggplot(plot_data, aes(x = Outcome, y = Parameter_Label, fill = original)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", original)), color = "black", size = 3) +
  scale_fill_gradient2(
    name = "PRCC",
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", # Blue-White-Red Diverging
    midpoint = 0, limits = c(-1, 1)
  ) +
  labs(title = "Global Sensitivity Analysis (LHS-PRCC)",
       x = "Outcomes (10-Year Cumulative)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(heatmap_plot)
ggsave("Figure_PRCC_Heatmap.png", heatmap_plot, width = 10, height = 12, dpi = 300)

