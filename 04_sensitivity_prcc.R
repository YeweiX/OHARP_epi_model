# =============================================================================
#
# Global Sensitivity Analysis using LHS-PRCC
#
# Description: This script performs a global sensitivity analysis on a 16-state
#              epidemiological model. It uses Latin Hypercube Sampling (LHS)
#              from the 'lhs' package to generate parameter sets, runs the model
#              in parallel, and then calculates Partial Rank Correlation
#              Coefficients (PRCC) to quantify parameter impact on key outcomes.
#
# =============================================================================


# --- 1. Load Necessary Libraries ---
# -----------------------------------------------------------------------------
library(deSolve)      # For solving Ordinary Differential Equations (ODEs)
library(lhs)          # For Latin Hypercube Sampling (LHS) -> Corrected
library(sensitivity)  # For PRCC (pcc) function
library(dplyr)        # For data manipulation and pipelines
library(ggplot2)      # For creating high-quality plots
library(doParallel)   # For enabling parallel processing


# --- 2. Setup Parallel Processing ---
# -----------------------------------------------------------------------------
# Use a reasonable number of cores, leaving some for system stability
num_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Parallel processing enabled on", num_cores, "cores.\n")


# --- 3. Define Model, Initial Conditions, and Time ---
# -----------------------------------------------------------------------------
# Set a seed for reproducibility of the random sampling process
set.seed(123)

# --- Initial Conditions (Latest Version) ---
N <- 5.18e6
Nh <- 37588 * 0.857
Nc <- N - Nh
RCh <- round(0.124 * Nh); RIh <- round(0.119 * 0.314 * 0.377 * Nh); SIh <- round(0.119 * 0.314 * 0.56 * Nh);
Sh <- Nh - RCh - RIh - SIh;
RCc <- round(0.257 * Nc); RIc <- round(0.001 * Nc); SIc <- round(0.01 * Nc);
Sc <- Nc - RCc - RIc - SIc;

# State vector for the ODE solver
y_init <- c(Sh = Sh, RCh = RCh, RIh = RIh, SIh = SIh, Sc = Sc, RCc = RCc, RIc = RIc, SIc = SIc,
            Drh = 0, Dsh = 0, Drc = 0, Dsc = 0, 
            C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0)

# --- Simulation Time ---
# We will analyze sensitivity at the end of a 10-year period
t_start <- 0
t_end <- 365 * 20
times <- seq(t_start, t_end, by = 365) # Annual output is sufficient
final_time <- max(times)

# --- ODE Model (Latest Version with Constrained Flows) ---
constrain_flow <- function(flow_rate, population, dt = 1) {
  max_flow <- population / dt
  return(min(flow_rate, max_flow))
}

infection_model_ode <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Ensure state variables are non-negative
    Sh <- max(0, y["Sh"]); RCh <- max(0, y["RCh"]); RIh <- max(0, y["RIh"]); SIh <- max(0, y["SIh"]);
    Sc <- max(0, y["Sc"]); RCc <- max(0, y["RCc"]); RIc <- max(0, y["RIc"]); SIc <- max(0, y["SIc"]);
    
    # Calculate current population sizes for density-dependent terms
    Nh_current <- max(1, Sh + RCh + RIh + SIh)
    Nc_current <- max(1, Sc + RCc + RIc + SIc)
    N_total_living <- Nh_current + Nc_current
    
    # Calculate all flows using the constrain_flow helper function for robustness
    flow_beta_eh_Sh <- constrain_flow(beta_eh * Sh, Sh); flow_beta_hh_Sh <- constrain_flow(beta_hh * Sh * RCh / Nh_current, Sh);
    flow_beta_sh_sih_Sh <- constrain_flow(beta_sh_sih * Sh * SIh / Nh_current, Sh); flow_alpha_dis_Sh <- constrain_flow(alpha_dis * Sh, Sh); flow_mu_d_Sh <- constrain_flow(mu_d * Sh, Sh);
    flow_delta_rch_rih_RCh <- constrain_flow(delta_rch_rih * RCh, RCh); flow_gamma_rch_sh_RCh <- constrain_flow(gamma_rch_sh * RCh, RCh);
    flow_alpha_dis_RCh <- constrain_flow(alpha_dis * RCh, RCh); flow_delta_rch_sih_RCh <- constrain_flow(delta_rch_sih * RCh, RCh); flow_mu_d_RCh <- constrain_flow(mu_d * RCh, RCh);
    flow_gamma_rih_rch_RIh <- constrain_flow(gamma_rih_rch * RIh, RIh); flow_mu_rih_RIh <- constrain_flow(mu_rih * RIh, RIh);
    flow_gamma_sih_sh_SIh <- constrain_flow(gamma_sih_sh * SIh, SIh); flow_mu_sih_SIh <- constrain_flow(mu_sih * SIh, SIh);
    flow_beta_sc_ec_Sc <- constrain_flow(beta_sc_ec * Sc, Sc); flow_beta_sc_ac_Sc <- constrain_flow(beta_sc_ac * Sc * 0.33, Sc);
    flow_beta_hc_Sc <- constrain_flow(beta_hc * Sc * RCc / Nc_current, Sc); flow_beta_sc_sic_Sc <- constrain_flow(beta_sc_sic * Sc * SIc / Nc_current, Sc);
    flow_alpha_adm_Sc <- constrain_flow(alpha_adm * Sc, Sc); flow_mu_d_Sc <- constrain_flow(mu_d * Sc, Sc);
    flow_delta_rcc_ric_RCc <- constrain_flow(delta_rcc_ric * RCc, RCc); flow_delta_rcc_sic_RCc <- constrain_flow(delta_rcc_sic * RCc, RCc);
    flow_gamma_rcc_sc_RCc <- constrain_flow(gamma_rcc_sc * RCc, RCc); flow_alpha_adm_RCc <- constrain_flow(alpha_adm * RCc, RCc); flow_mu_d_RCc <- constrain_flow(mu_d * RCc, RCc);
    flow_gamma_ric_rcc_RIc <- constrain_flow(gamma_ric_rcc * RIc, RIc); flow_mu_ric_RIc <- constrain_flow(mu_ric * RIc, RIc); flow_alpha_adm_RIc <- constrain_flow(alpha_adm * RIc, RIc);
    flow_gamma_sic_sc_SIc <- constrain_flow(gamma_sic_sc * SIc, SIc); flow_mu_sic_SIc <- constrain_flow(mu_sic * SIc, SIc); flow_alpha_adm_SIc <- constrain_flow(alpha_adm * SIc, SIc);
    
    # Define the differential equations
    dSh_dt <- -flow_beta_eh_Sh - flow_beta_hh_Sh - flow_beta_sh_sih_Sh - flow_alpha_dis_Sh - flow_mu_d_Sh + flow_gamma_sih_sh_SIh + flow_gamma_rch_sh_RCh + flow_alpha_adm_Sc;
    dRCh_dt <- -flow_delta_rch_rih_RCh - flow_gamma_rch_sh_RCh - flow_alpha_dis_RCh - flow_delta_rch_sih_RCh - flow_mu_d_RCh + flow_beta_eh_Sh + flow_beta_hh_Sh + flow_gamma_rih_rch_RIh + flow_alpha_adm_RCc;
    dRIh_dt <- -flow_gamma_rih_rch_RIh - flow_mu_rih_RIh + flow_delta_rch_rih_RCh + flow_alpha_adm_RIc;
    dSIh_dt <- -flow_gamma_sih_sh_SIh - flow_mu_sih_SIh + flow_beta_sh_sih_Sh + flow_delta_rch_sih_RCh + flow_alpha_adm_SIc;
    dSc_dt <- -flow_beta_sc_ec_Sc - flow_beta_sc_ac_Sc - flow_beta_hc_Sc - flow_beta_sc_sic_Sc - flow_alpha_adm_Sc - flow_mu_d_Sc + flow_gamma_sic_sc_SIc + flow_gamma_rcc_sc_RCc + flow_alpha_dis_Sh + mu_b * N_total_living + mu_m * N_total_living;
    dRCc_dt <- -flow_delta_rcc_ric_RCc - flow_delta_rcc_sic_RCc - flow_gamma_rcc_sc_RCc - flow_alpha_adm_RCc - flow_mu_d_RCc + flow_beta_sc_ec_Sc + flow_beta_sc_ac_Sc + flow_beta_hc_Sc + flow_gamma_ric_rcc_RIc + flow_alpha_dis_RCh;
    dRIc_dt <- -flow_gamma_ric_rcc_RIc - flow_mu_ric_RIc - flow_alpha_adm_RIc + flow_delta_rcc_ric_RCc;
    dSIc_dt <- -flow_gamma_sic_sc_SIc - flow_mu_sic_SIc - flow_alpha_adm_SIc + flow_beta_sc_sic_Sc + flow_delta_rcc_sic_RCc;
    dDrh_dt <- flow_mu_rih_RIh; dDsh_dt <- flow_mu_sih_SIh; dDrc_dt <- flow_mu_ric_RIc; dDsc_dt <- flow_mu_sic_SIc;
    dC_RIh_dt <- flow_delta_rch_rih_RCh + flow_alpha_adm_RIc; dC_SIh_dt <- flow_beta_sh_sih_Sh + flow_delta_rch_sih_RCh + flow_alpha_adm_SIc;
    dC_RIc_dt <- delta_rcc_ric * RCc; dC_SIc_dt <- flow_beta_sc_sic_Sc + flow_delta_rcc_sic_RCc;
    
    # Return the derivatives as a list
    return(list(c(dSh_dt, dRCh_dt, dRIh_dt, dSIh_dt, dSc_dt, dRCc_dt, dRIc_dt, dSIc_dt, dDrh_dt, dDsh_dt, dDrc_dt, dDsc_dt, dC_RIh_dt, dC_SIh_dt, dC_RIc_dt, dC_SIc_dt)))
  })
}


# --- 4. Define Parameter Ranges and Generate Samples with LHS ---
# -----------------------------------------------------------------------------
cat("\nGenerating parameter samples using Latin Hypercube Sampling...\n")


#parameter distributions

epi_param_distributions <- list(
  beta_hh       = list(dist = "unif", min = 0.002100771, max = 0.013899509),
  beta_eh       = list(dist = "unif", min =0.000000084,  max = 0.033214594),
  beta_hc       = list(dist = "unif", min = 0.001600901, max = 0.010231657),
  beta_sc_ac    = list(dist = "unif", min = 0.000000096, max = 0.003083009),
  beta_sc_ec    = list(dist = "unif", min = 0.000000037, max = 0.001111809),  # was 0 → use eps
  beta_sh_sih   = list(dist = "unif", min = 0.155454327, max = 0.996909995),
  delta_rcc_ric = list(dist = "unif", min =0.000031951,  max = 0.000149248),
  
  # --- For the other parameters, we use the ranges you originally defined ---
  delta_rch_rih = list(dist = "unif", min = 0.00010513, max = 0.06779661),
  mu_sih        = list(dist = "unif", min = 0.00233333, max = 0.00390000),
  mu_rih        = list(dist = "unif", min = 0.00254333, max = 0.00554667),
  gamma_sih_sh  = list(dist = "unif", min = 0.09090909, max = 0.20000000),
  gamma_rch_sh  = list(dist = "unif", min = 0.00147397, max = 0.01076667),
  gamma_rih_rch = list(dist = "unif", min = 0.06666667, max = 0.16666667), #1/los
  delta_rch_sih = list(dist = "unif", min = 0.000592603,max = 0.001100548),
  beta_sc_sic   = list(dist = "unif", min = 0.00000707, max = 0.00000877),
  gamma_sic_sc  = list(dist = "unif", min = 0.16666667, max = 0.33333333),
  gamma_rcc_sc  = list(dist = "unif", min = 0.00460000, max = 0.01800000),
  gamma_ric_rcc = list(dist = "unif", min = 0.12500000, max = 0.25000000),
  mu_sic        = list(dist = "unif", min = 0.003226*0.7, max = 0.003226*1.3),
  mu_ric        = list(dist = "unif", min = 0.003629*0.7, max = 0.003629*1.3),
  delta_rcc_sic = list(dist = "unif", min = 0.000592603,max = 0.001100548),
  alpha_adm     = list(dist = "unif", min = 0.000316 * 0.7, max = 0.000316 * 1.3),
  alpha_dis     = list(dist = "unif", min = 0.235 * 0.7,    max = 0.235 * 1.3),
  mu_b          = list(dist = "unif", min = 0.00002518 * 0.7, max = 0.00002518 * 1.3),
  mu_m          = list(dist = "unif", min = 0.000002 * 0.7,   max = 0.000002 * 1.3),
  mu_d          = list(dist = "unif", min = 0.00001353 * 0.7, max = 0.00001353 * 1.3)
)



# Extract parameter names and set number of samples
param_names_to_vary <- names(epi_param_distributions)
n_params <- length(param_names_to_vary)
n_samples <- 1000 # Number of simulations (increase for higher precision, e.g., 1000)

# --- CORRECTED: Use randomLHS from the 'lhs' package ---
# This generates a matrix of probabilities (0-1) for each parameter
lhs_matrix <- randomLHS(n = n_samples, k = n_params)

# Create an empty data frame for parameter values
param_samples <- data.frame(matrix(NA, nrow = n_samples, ncol = n_params))
colnames(param_samples) <- param_names_to_vary

# Convert LHS probabilities to parameter values using their respective distributions
for (i in 1:n_params) {
  param_name <- param_names_to_vary[i]
  dist_info <- epi_param_distributions[[param_name]]
  
  if (dist_info$dist == "unif") {
    # Use the quantile function for the uniform distribution (qunif)
    param_samples[, i] <- qunif(lhs_matrix[, i], min = dist_info$min, max = dist_info$max)
  }
  # Note: Add 'else if' blocks here for other distributions (e.g., qnorm, qbeta) if needed.
}
summary(param_samples)



economic_horizon_years <- 10
full_sim_years <- 20 # Assuming a 20-year total simulation
economic_start_year <- full_sim_years - economic_horizon_years
# --- 5. Run Simulations in Parallel ---
# -----------------------------------------------------------------------------
cat("\nRunning", n_samples, "simulations in parallel...\n")

results_list <- foreach(i = 1:n_samples, .packages = c("deSolve", "dplyr")) %dopar% {
  current_params <- as.list(param_samples[i, ])
  
  # Run the ODE solver
  out <- try(ode(y      = y_init,
                 times  = times,
                 func   = infection_model_ode,
                 parms  = current_params,
                 method = "lsoda"),
             silent = TRUE)
  if (inherits(out, "try-error")) {
    return(NULL)
  }
  
  out_df <- as.data.frame(out)
  
  # Define start of the last-10-year window
  start_time <- economic_start_year * 365
  start_idx  <- which(out_df$time >= start_time)[1]
  if (is.na(start_idx) || start_idx == nrow(out_df)) {
    return(NULL)
  }
  
  # Slice to the last 10 years and compute deltas
  recent      <- out_df[start_idx:nrow(out_df), ]
  delta_C_RIh <- tail(recent$C_RIh, 1) - recent$C_RIh[1]
  delta_C_RIc <- tail(recent$C_RIc, 1) - recent$C_RIc[1]
  delta_Drh   <- tail(recent$Drh,   1) - recent$Drh[1]
  delta_Drc   <- tail(recent$Drc,   1) - recent$Drc[1]
  
  data.frame(
    sample_id = i,
    C_RIh     = delta_C_RIh,
    C_RIc     = delta_C_RIc,
    Drh  = delta_Drh,
    Drc  = delta_Drc
  )
}


  

# Stop the parallel cluster
stopCluster(cl)
cat("Simulations complete.\n")

# Combine results 
all_results <- bind_rows(results_list)
cat(nrow(all_results), "out of", n_samples, "simulations completed successfully.\n")

# --- 6. Calculate PRCC for Key Outcomes ---
# -----------------------------------------------------------------------------
cat("\nCalculating PRCC for key epidemiological outcomes...\n")

outcome_vars <- c("C_RIh", "C_RIc", "Drh", "Drc")

outcome_titles <- c(
  C_RIh = "Cumulative Hospital Resistant Infection",
  C_RIc = "Cumulative Community Resistant Infection",
  Drh   = "Deaths from Hospital Resistant Infection",
  Drc   = "Deaths from Community Resistant Infection"
)

# Keep only successful runs
X <- param_samples[all_results$sample_id, , drop = FALSE]

# Safety check
if (nrow(X) != nrow(all_results)) {
  stop("Mismatch between parameter samples and simulation results.")
}

# Function to extract PRCC results from sensitivity::pcc()
extract_prcc_df <- function(prcc_result, outcome_name) {
  
  prcc_mat <- as.data.frame(prcc_result$PRCC)
  
  prcc_mat$Parameter <- rownames(prcc_mat)
  rownames(prcc_mat) <- NULL
  
  # Standardise column names in case R changes punctuation
  names(prcc_mat) <- gsub(" ", "_", names(prcc_mat))
  names(prcc_mat) <- gsub("\\.", "_", names(prcc_mat))
  
  # Rename original PRCC column if needed
  if (!"original" %in% names(prcc_mat)) {
    original_col <- names(prcc_mat)[sapply(prcc_mat, is.numeric)][1]
    names(prcc_mat)[names(prcc_mat) == original_col] <- "original"
  }
  
  prcc_mat$Outcome <- outcome_name
  
  return(prcc_mat)
}

prcc_results_list <- list()
all_prcc_df <- list()

for (outcome_var in outcome_vars) {
  
  outcome_name <- outcome_titles[[outcome_var]]
  y <- as.numeric(all_results[[outcome_var]])
  
  # Remove non-finite values
  valid_idx <- is.finite(y) & complete.cases(X)
  X_valid <- X[valid_idx, , drop = FALSE]
  y_valid <- y[valid_idx]
  
  # Avoid PRCC failure if outcome is constant
  if (length(unique(y_valid)) <= 1) {
    warning(paste("Skipping", outcome_name, "because outcome is constant."))
    next
  }
  
  prcc_result <- pcc(
    X = X_valid,
    y = y_valid,
    rank = TRUE,
    nboot = 1000
  )
  
  prcc_results_list[[outcome_name]] <- prcc_result
  all_prcc_df[[outcome_name]] <- extract_prcc_df(prcc_result, outcome_name)
}

all_prcc_df <- bind_rows(all_prcc_df)

cat("PRCC calculations finished.\n")


# --- 7. Prepare PRCC Results for Plotting ---
# -----------------------------------------------------------------------------

param_full_names <- c(
  # Hospital Acquisition & Progression
  "beta_hh"       = "Acquisition Rate (Hospital): Human Contact",
  "beta_sh_sih"   = "Acquisition Rate (Hospital): Susceptible to Infection",
  "beta_eh"       = "Acquisition Rate (Hospital): Environment",
  "delta_rch_rih" = "Progression Rate (Hospital): Colonised to Resistant Infection",
  "delta_rch_sih" = "Progression Rate (Hospital): Colonised to Susceptible Infection",
  
  # Hospital Recovery & Mortality
  "gamma_sih_sh"  = "Recovery Rate (Hospital): Infection to Susceptible",
  "gamma_rih_rch" = "Recovery Rate (Hospital): Resistant Infection to Colonised",
  "gamma_rch_sh"  = "Decolonisation Rate (Hospital): Colonised to Susceptible",
  "mu_sih"        = "Mortality Rate (Hospital): Susceptible Infection",
  "mu_rih"        = "Mortality Rate (Hospital): Resistant Infection",
  
  # Community Acquisition & Progression
  "beta_hc"       = "Acquisition Rate (Community): Human Contact",
  "beta_sc_sic"   = "Acquisition Rate (Community): Susceptible to Infection",
  "beta_sc_ac"    = "Acquisition Rate (Community): Animal Contact",
  "beta_sc_ec"    = "Acquisition Rate (Community): Environment",
  "delta_rcc_ric" = "Progression Rate (Community): Colonised to Resistant Infection",
  "delta_rcc_sic" = "Progression Rate (Community): Colonised to Susceptible Infection",
  
  # Community Recovery & Mortality
  "gamma_sic_sc"  = "Recovery Rate (Community): Infection to Susceptible",
  "gamma_ric_rcc" = "Recovery Rate (Community): Resistant Infection to Colonised",
  "gamma_rcc_sc"  = "Decolonisation Rate (Community): Colonised to Susceptible",
  "mu_sic"        = "Mortality Rate (Community): Susceptible Infection",
  "mu_ric"        = "Mortality Rate (Community): Resistant Infection",
  
  # General Demographics & Hospital Flow
  "alpha_adm"     = "Hospital Admission Rate",
  "alpha_dis"     = "Hospital Discharge Rate",
  "mu_b"          = "Crude Birth Rate",
  "mu_m"          = "Net Migration Rate",
  "mu_d"          = "Crude Death Rate"
)

prcc_data_for_plot <- all_prcc_df %>%
  mutate(
    Parameter_Full = param_full_names[Parameter],
    Outcome = factor(Outcome, levels = unname(outcome_titles))
  ) %>%
  filter(!is.na(Parameter_Full)) %>%
  filter(is.finite(original))

write.csv(prcc_data_for_plot, "prcc_results.csv", row.names = FALSE)


# --- 8. Generate PRCC Heatmap ---
# -----------------------------------------------------------------------------

cat("Generating PRCC heatmap...\n")

prcc_heatmap_final <- ggplot(
  prcc_data_for_plot,
  aes(x = Outcome, y = Parameter_Full, fill = original)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", original)),
    color = "white",
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    name = "PRCC",
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(
    title = "PRCC Global Sensitivity Analysis",
    x = "Model Outcome",
    y = "Input Parameter"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

print(prcc_heatmap_final)

ggsave(
  filename = "prcc_heatmap.png",
  plot = prcc_heatmap_final,
  width = 11,
  height = 8,
  dpi = 300
)


# --- 9. Optional Monotonicity Check ---
# -----------------------------------------------------------------------------

check_monotonicity <- function(X, y, outcome_name, n_top = 6, plot = TRUE) {
  
  cor_vec <- apply(X, 2, function(col) {
    cor(col, y, method = "spearman", use = "complete.obs")
  })
  
  top_params <- names(sort(abs(cor_vec), decreasing = TRUE))[1:n_top]
  
  if (plot) {
    plot_list <- lapply(top_params, function(param) {
      ggplot(data.frame(x = X[[param]], y = y), aes(x = x, y = y)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "loess", se = FALSE) +
        labs(
          title = paste("Monotonicity Check:", param),
          subtitle = paste("Outcome:", outcome_name),
          x = param,
          y = outcome_name
        ) +
        theme_minimal()
    })
    
    n_col <- min(n_top, 3)
    do.call(gridExtra::grid.arrange, c(plot_list, ncol = n_col))
  }
  
  return(cor_vec[top_params])
}

# Example monotonicity check for the final outcome
last_outcome_var <- "Drc"
last_outcome_name <- outcome_titles[[last_outcome_var]]
last_y <- as.numeric(all_results[[last_outcome_var]])

cat("\nMonotonicity check for outcome:", last_outcome_name, "\n")
top_correlations <- check_monotonicity(
  X = X,
  y = last_y,
  outcome_name = last_outcome_name,
  n_top = 6,
  plot = TRUE
)

print(round(top_correlations, 3))

