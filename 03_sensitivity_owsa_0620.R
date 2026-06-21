# =============================================================================
#
# One-Way Sensitivity Analysis (OWSA) for an Epidemiological Model
#
# =============================================================================


# --- 1. SETUP & PREREQUISITES ---
# -----------------------------------------------------------------------------
# Ensure all required packages are loaded ONCE at the top.
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)

# --- Analysis Settings ---
# NOTE: The following objects must be defined in your R environment before running:
# 1. infection_model_ode: The function defining your ODE model.
# 2. y_init_baseline: A named vector of the model's initial conditions.
# 3. param_full_names: A named vector mapping short parameter names to full names for plots.
#
# Example: source("your_model_setup_file.R")

economic_horizon_years <- 10
full_sim_years <- 20
economic_start_year <- full_sim_years - economic_horizon_years

# Define the baseline parameter values
baseline_epi_params <- list(
  beta_hh =0.007893581, beta_sh_sih = 0.239527636, delta_rch_rih = 0.03395087,
  mu_sih = 0.00311667, mu_rih = 0.00404500, gamma_sih_sh = 0.14285714,
  gamma_rch_sh = 0.00612032, gamma_rih_rch = 0.11111111, delta_rch_sih = 0.00084658,
  beta_sc_sic = 0.00000792, delta_rcc_ric = 0.000085556, delta_rcc_sic = 0.00084658,
  gamma_sic_sc = 0.25, gamma_rcc_sc = 0.009, gamma_ric_rcc = 0.17,
  mu_sic =  0.003226, mu_ric =  0.003629, alpha_adm = 0.000316, alpha_dis = 0.235,
  beta_eh = 0.006259864, beta_hc = 0.003655865, beta_sc_ac =  0.000683476, beta_sc_ec = 0.000260089,
  mu_b = 0.00002518, mu_m = 0.000002, mu_d = 0.00001353
)


# --- 2. HELPER FUNCTION TO RUN SIMULATIONS ---
# -----------------------------------------------------------------------------
# This function runs one simulation and correctly extracts the total outcomes.
run_and_extract_outcomes <- function(parms) {
  times <- seq(0, full_sim_years * 365, by = 365)
  start_time <- economic_start_year * 365
  
  # Run the ODE solver
  out <- try(ode(y = y_init, times = times, func = infection_model_ode, parms = parms, method = "lsoda"), silent = TRUE)
  
  # Handle solver errors
  if (inherits(out, "try-error")) {
    return(list(Drh=NA, Drc=NA, Cases_h=NA, Cases_c=NA))
  }
  
  out_df <- as.data.frame(out)
  start_idx <- which(out_df$time >= start_time)[1]
  
  # Handle cases where simulation did not run long enough
  if (is.na(start_idx) || start_idx >= nrow(out_df)) {
    return(list(Drh=NA, Drc=NA, Cases_h=NA, Cases_c=NA))
  }
  
  # Get the values at the start and end of the economic horizon
  first_row <- out_df[start_idx, ]
  last_row <- out_df[nrow(out_df), ]
  
  # Calculate the change (delta) for each outcome
  total_Drh <- last_row$Drh - first_row$Drh
  total_Drc <- last_row$Drc - first_row$Drc
  total_Cases_h <- (last_row$C_RIh) - (first_row$C_RIh)
  total_Cases_c <- (last_row$C_RIc) - (first_row$C_RIc)
  
  return(list(Drh = total_Drh, Drc = total_Drc, Cases_h = total_Cases_h, Cases_c = total_Cases_c))
}


# --- 3. PRE-ANALYSIS: IDENTIFY TOP PARAMETERS ---
# -----------------------------------------------------------------------------
cat("Step 1: Identifying influential parameters via pre-scan...\n")

excluded_params <- c("alpha_adm", "alpha_dis")

owsa_pre_scan_list <- list()
for (param_name in names(baseline_epi_params)) {
  if (param_name %in% excluded_params) next
  
  baseline_val <- baseline_epi_params[[param_name]]
  if (is.na(baseline_val) || baseline_val == 0) next
  
  params_low <- baseline_epi_params; params_low[[param_name]] <- baseline_val * 0.5
  out_low <- run_and_extract_outcomes(params_low)
  
  params_high <- baseline_epi_params; params_high[[param_name]] <- baseline_val * 1.5
  out_high <- run_and_extract_outcomes(params_high)
  
  owsa_pre_scan_list[[param_name]] <- data.frame(
    parameter = param_name,
    range_Drh = abs(out_high$Drh - out_low$Drh),
    range_Drc = abs(out_high$Drc - out_low$Drc),
    range_Cases_h = abs(out_high$Cases_h - out_low$Cases_h),
    range_Cases_c = abs(out_high$Cases_c - out_low$Cases_c)
  )
}
owsa_pre_scan_df <- bind_rows(owsa_pre_scan_list) %>% filter(if_any(starts_with("range"), ~ !is.na(.)))

# Identify the Top 5 most influential parameters for each outcome
top5_Drh     <- owsa_pre_scan_df %>% arrange(desc(range_Drh))     %>% slice_head(n = 5) %>% pull(parameter)
top5_Drc     <- owsa_pre_scan_df %>% arrange(desc(range_Drc))     %>% slice_head(n = 5) %>% pull(parameter)
top5_Cases_h <- owsa_pre_scan_df %>% arrange(desc(range_Cases_h)) %>% slice_head(n = 5) %>% pull(parameter)
top5_Cases_c <- owsa_pre_scan_df %>% arrange(desc(range_Cases_c)) %>% slice_head(n = 5) %>% pull(parameter)

# Create a final, unique list of parameters for the main simulation
key_parameters <- unique(c(top5_Drh, top5_Drc, top5_Cases_h, top5_Cases_c))
cat("Pre-scan complete. Key parameters for detailed analysis:", paste(key_parameters, collapse=", "), "\n")


# --- 4. MAIN SIMULATION: PARAMETER SWEEP ---
# -----------------------------------------------------------------------------
cat("Step 2: Running main simulation sweep...\n")

percent_changes <- seq(-1, 1, by = 0.25)
baseline_outcomes <- run_and_extract_outcomes(baseline_epi_params)
heatmap_data_list <- list()

for (param_name in key_parameters) {
  for (pc in percent_changes) {
    if (pc == 0) next
    current_params <- baseline_epi_params
    current_params[[param_name]] <- current_params[[param_name]] * (1 + pc)
    outcomes <- run_and_extract_outcomes(current_params)
    
    heatmap_data_list[[length(heatmap_data_list) + 1]] <- data.frame(
      parameter = param_name,
      percent_change = pc,
      change_Drh = outcomes$Drh - baseline_outcomes$Drh,
      change_Drc = outcomes$Drc - baseline_outcomes$Drc,
      change_Cases_h = outcomes$Cases_h - baseline_outcomes$Cases_h,
      change_Cases_c = outcomes$Cases_c - baseline_outcomes$Cases_c
    )
  }
}
heatmap_df <- bind_rows(heatmap_data_list)
cat("Simulations complete.\n")


# --- 5. DATA PREPARATION AND PLOTTING ---
# -----------------------------------------------------------------------------
cat("Step 3: Preparing data and generating plots...\n")

# Pivot the data to a long format suitable for ggplot
heatmap_long <- heatmap_df %>%
  pivot_longer(
    cols = starts_with("change_"),
    names_to = "outcome_name",
    values_to = "absolute_change"
  ) %>%
  mutate(
    type   = ifelse(grepl("^change_Dr", outcome_name), "Deaths", "Cases"),
    origin = ifelse(grepl("h$", outcome_name), "Hospital", "Community")
  )

# Filter the data for the two plots, ensuring only top parameters are shown for each
deaths_plot_data_filtered <- heatmap_long %>%
  filter(
    type == "Deaths",
    (origin == "Hospital"  & parameter %in% top5_Drh) |
      (origin == "Community" & parameter %in% top5_Drc)
  )

cases_plot_data_filtered <- heatmap_long %>%
  filter(
    type == "Cases",
    (origin == "Hospital"  & parameter %in% top5_Cases_h) |
      (origin == "Community" & parameter %in% top5_Cases_c)
  )


