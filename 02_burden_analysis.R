# =============================================================================
# Scenario Analysis: Baseline vs. No Resistant Infections
# Calculating incremental burden (Baseline - No Resistance Scenario)
# Outcomes: Annual Drh, Drc deaths and cumulative RIh, RIc cases
# =============================================================================
load("/mcmc.RData") ##load calibration result here

# --- 0. Load Necessary Libraries ---
library(deSolve)
library(dplyr)
library(lhs)
library(sensitivity)
library(purrr)
library(ggplot2)
library(scales)
library(rstan)
library(doParallel)
library(tidyr)

# --- Setup Parallel Processing ---
num_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Parallel processing enabled on", num_cores, "cores.\n")

# --- 1. Load MCMC Posterior Samples ---
mcmc_param_names <- c("beta_hh", "beta_eh", "beta_hc", "beta_sc_ac", "beta_sc_ec", "beta_sh_sih", "delta_rcc_ric")

# Baseline fixed parameters
params <- list(
  beta_hh =0.007893581, beta_sh_sih = 0.239527636, delta_rch_rih = 0.03395087,
  mu_sih = 0.00311667, mu_rih = 0.00404500, gamma_sih_sh = 0.14285714,
  gamma_rch_sh = 0.00612032, gamma_rih_rch = 0.11111111, delta_rch_sih = 0.00084658,
  beta_sc_sic = 0.00000792, delta_rcc_ric = 0.000085556, delta_rcc_sic = 0.00084658,
  gamma_sic_sc = 0.25, gamma_rcc_sc = 0.009, gamma_ric_rcc = 0.17,
  mu_sic =  0.003226, mu_ric =  0.003629, alpha_adm = 0.000316, alpha_dis = 0.235,
  beta_eh = 0.006259864, beta_hc = 0.003655865, beta_sc_ac =  0.000683476, beta_sc_ec = 0.000260089,
  mu_b = 0.00002518, mu_m = 0.000002, mu_d = 0.00001353
)

baseline_fixed_params <- all_model_params_list[!names(all_model_params_list) %in% mcmc_param_names]

# Extract MCMC posterior samples
posterior_samples_list <- rstan::extract("xxx", pars = mcmc_param_names, permuted = TRUE) ##use calibration result xxx here

stan_list_to_df <- function(stan_list) {
  dfs <- list()
  for (name in names(stan_list)) {
    samples <- stan_list[[name]]
    if (is.matrix(samples)) {
      temp_df <- as.data.frame(t(samples))
      colnames(temp_df) <- paste0(name, "[", 1:ncol(temp_df), "]")
      dfs[[name]] <- temp_df
    } else {
      dfs[[name]] <- data.frame(samples)
      colnames(dfs[[name]]) <- name
    }
  }
  final_df <- bind_cols(dfs)
  return(final_df)
}

posterior_samples_mcmc_selected <- stan_list_to_df(posterior_samples_list)
cat("Extracted MCMC parameters successfully.\n")
cat("Shape of MCMC samples dataframe:", dim(posterior_samples_mcmc_selected), "\n")

# --- 2. Define Initial Conditions ---
N <- 5.18e6
Nh <- 37588 * 0.857
Nc <- N - Nh

# BASELINE Initial Conditions (with existing resistant infections)
RCh_baseline <- round(0.124 * Nh)
RIh_baseline <- round(0.119 * 0.314 * 0.377 * Nh)
SIh_baseline <- round(0.119 * 0.314 * 0.56 * Nh)
Sh_baseline <- Nh - RCh_baseline - RIh_baseline - SIh_baseline

RCc_baseline <- round(0.257 * Nc)
RIc_baseline <- round(0.001 * Nc)
SIc_baseline <- round(0.01 * Nc)
Sc_baseline <- Nc - RCc_baseline - RIc_baseline - SIc_baseline

y_init_baseline <- c(
  Sh = Sh_baseline, RCh = RCh_baseline, RIh = RIh_baseline, SIh = SIh_baseline,
  Sc = Sc_baseline, RCc = RCc_baseline, RIc = RIc_baseline, SIc = SIc_baseline,
  Drh = 0, Dsh = 0, Drc = 0, Dsc = 0,
  C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0
)

# NO-RESISTANCE Initial Conditions (RIh = 0, RIc = 0)
# Redistribute resistant infections back to susceptible compartments
RCh_no_res <- RCh_baseline
RIh_no_res <- 0  # No resistant infections
SIh_no_res <- SIh_baseline
Sh_no_res <- Nh - RCh_no_res - RIh_no_res - SIh_no_res  # Absorb the RIh into Sh

RCc_no_res <- RCc_baseline
RIc_no_res <- 0  # No resistant infections
SIc_no_res <- SIc_baseline
Sc_no_res <- Nc - RCc_no_res - RIc_no_res - SIc_no_res  # Absorb the RIc into Sc

y_init_no_resistance <- c(
  Sh = Sh_no_res, RCh = RCh_no_res, RIh = RIh_no_res, SIh = SIh_no_res,
  Sc = Sc_no_res, RCc = RCc_no_res, RIc = RIc_no_res, SIc = SIc_no_res,
  Drh = 0, Dsh = 0, Drc = 0, Dsc = 0,
  C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0
)

cat("Initial conditions set up:\n")
cat("Baseline - RIh:", RIh_baseline, "RIc:", RIc_baseline, "\n")
cat("No-resistance - RIh:", RIh_no_res, "RIc:", RIc_no_res, "\n")

# --- 3. Time settings ---
full_sim_years <- 20
economic_horizon_years <- 10

t_start <- 0
t_end <- 365 * full_sim_years
times <- seq(t_start, t_end, by = 365)  # Annual time points

# --- 4. ODE Model (same as before) ---
constrain_flow <- function(flow_rate, population, dt = 1) {
  max_flow <- population / dt
  return(min(flow_rate, max_flow))
}

infection_model_ode <- function(t, y, params) {
  with(as.list(c(y, params)), {
    Sh <- max(0, y["Sh"]); RCh <- max(0, y["RCh"]); RIh <- max(0, y["RIh"]); SIh <- max(0, y["SIh"]);
    Sc <- max(0, y["Sc"]); RCc <- max(0, y["RCc"]); RIc <- max(0, y["RIc"]); SIc <- max(0, y["SIc"]);
    
    Nh_current <- max(1, Sh + RCh + RIh + SIh)
    Nc_current <- max(1, Sc + RCc + RIc + SIc)
    N_total_living <- Nh_current + Nc_current
    
    # Hospital flows OUT
    flow_beta_eh_Sh <- constrain_flow(beta_eh * Sh, Sh)
    flow_beta_hh_Sh <- constrain_flow(beta_hh * Sh * RCh / Nh_current, Sh)
    flow_beta_sh_sih_Sh <- constrain_flow(beta_sh_sih * Sh * SIh / Nh_current, Sh)
    flow_alpha_dis_Sh <- constrain_flow(alpha_dis * Sh, Sh)
    flow_mu_d_Sh <- constrain_flow(mu_d * Sh, Sh)
    
    flow_delta_rch_rih_RCh <- constrain_flow(delta_rch_rih * RCh, RCh)
    flow_gamma_rch_sh_RCh <- constrain_flow(gamma_rch_sh * RCh, RCh)
    flow_alpha_dis_RCh <- constrain_flow(alpha_dis * RCh, RCh)
    flow_delta_rch_sih_RCh <- constrain_flow(delta_rch_sih * RCh, RCh)
    flow_mu_d_RCh <- constrain_flow(mu_d * RCh, RCh)
    
    flow_gamma_rih_rch_RIh <- constrain_flow(gamma_rih_rch * RIh, RIh)
    flow_mu_rih_RIh <- constrain_flow(mu_rih * RIh, RIh)
    
    flow_gamma_sih_sh_SIh <- constrain_flow(gamma_sih_sh * SIh, SIh)
    flow_mu_sih_SIh <- constrain_flow(mu_sih * SIh, SIh)
    
    # Community flows OUT
    flow_beta_sc_ec_Sc <- constrain_flow(beta_sc_ec * Sc, Sc)
    flow_beta_sc_ac_Sc <- constrain_flow(beta_sc_ac * Sc * 0.33, Sc)
    flow_beta_hc_Sc <- constrain_flow(beta_hc * Sc * RCc / Nc_current, Sc)
    flow_beta_sc_sic_Sc <- constrain_flow(beta_sc_sic * Sc * SIc / Nc_current, Sc)
    flow_alpha_adm_Sc <- constrain_flow(alpha_adm * Sc, Sc)
    flow_mu_d_Sc <- constrain_flow(mu_d * Sc, Sc)
    
    flow_delta_rcc_ric_RCc <- constrain_flow(delta_rcc_ric * RCc, RCc)
    flow_delta_rcc_sic_RCc <- constrain_flow(delta_rcc_sic * RCc, RCc)
    flow_gamma_rcc_sc_RCc <- constrain_flow(gamma_rcc_sc * RCc, RCc)
    flow_alpha_adm_RCc <- constrain_flow(alpha_adm * RCc, RCc)
    flow_mu_d_RCc <- constrain_flow(mu_d * RCc, RCc)
    
    flow_gamma_ric_rcc_RIc <- constrain_flow(gamma_ric_rcc * RIc, RIc)
    flow_mu_ric_RIc <- constrain_flow(mu_ric * RIc, RIc)
    flow_alpha_adm_RIc <- constrain_flow(alpha_adm * RIc, RIc)
    
    flow_gamma_sic_sc_SIc <- constrain_flow(gamma_sic_sc * SIc, SIc)
    flow_mu_sic_SIc <- constrain_flow(mu_sic * SIc, SIc)
    flow_alpha_adm_SIc <- constrain_flow(alpha_adm * SIc, SIc)
    
    # Hospital compartments
    dSh_dt <- -flow_beta_eh_Sh - flow_beta_hh_Sh - flow_beta_sh_sih_Sh - flow_alpha_dis_Sh - flow_mu_d_Sh +
      flow_gamma_sih_sh_SIh + flow_gamma_rch_sh_RCh + flow_alpha_adm_Sc
    dRCh_dt <- -flow_delta_rch_rih_RCh - flow_gamma_rch_sh_RCh - flow_alpha_dis_RCh - flow_delta_rch_sih_RCh - flow_mu_d_RCh +
      flow_beta_eh_Sh + flow_beta_hh_Sh + flow_gamma_rih_rch_RIh + flow_alpha_adm_RCc
    dRIh_dt <- -flow_gamma_rih_rch_RIh - flow_mu_rih_RIh + flow_delta_rch_rih_RCh + flow_alpha_adm_RIc
    dSIh_dt <- -flow_gamma_sih_sh_SIh - flow_mu_sih_SIh + flow_beta_sh_sih_Sh + flow_delta_rch_sih_RCh + flow_alpha_adm_SIc
    
    # Community compartments
    dSc_dt <- -flow_beta_sc_ec_Sc - flow_beta_sc_ac_Sc - flow_beta_hc_Sc - flow_beta_sc_sic_Sc - flow_alpha_adm_Sc - flow_mu_d_Sc +
      flow_gamma_sic_sc_SIc + flow_gamma_rcc_sc_RCc + flow_alpha_dis_Sh + mu_b * N_total_living + mu_m * N_total_living
    dRCc_dt <- -flow_delta_rcc_ric_RCc - flow_delta_rcc_sic_RCc - flow_gamma_rcc_sc_RCc - flow_alpha_adm_RCc - flow_mu_d_RCc +
      flow_beta_sc_ec_Sc + flow_beta_sc_ac_Sc + flow_beta_hc_Sc + flow_gamma_ric_rcc_RIc + flow_alpha_dis_RCh
    dRIc_dt <- -flow_gamma_ric_rcc_RIc - flow_mu_ric_RIc - flow_alpha_adm_RIc + flow_delta_rcc_ric_RCc
    dSIc_dt <- -flow_gamma_sic_sc_SIc - flow_mu_sic_SIc - flow_alpha_adm_SIc + flow_beta_sc_sic_Sc + flow_delta_rcc_sic_RCc
    
    # Cumulative deaths
    dDrh_dt <- flow_mu_rih_RIh
    dDsh_dt <- flow_mu_sih_SIh
    dDrc_dt <- flow_mu_ric_RIc
    dDsc_dt <- flow_mu_sic_SIc
    
    # Cumulative new cases
    dC_RIh_dt <- flow_delta_rch_rih_RCh + flow_alpha_adm_RIc
    dC_SIh_dt <- flow_beta_sh_sih_Sh + flow_delta_rch_sih_RCh + flow_alpha_adm_SIc
    dC_RIc_dt <- flow_delta_rcc_ric_RCc
    dC_SIc_dt <- flow_beta_sc_sic_Sc + flow_delta_rcc_sic_RCc
    
    return(list(c(dSh_dt, dRCh_dt, dRIh_dt, dSIh_dt, dSc_dt, dRCc_dt, dRIc_dt, dSIc_dt,
                  dDrh_dt, dDsh_dt, dDrc_dt, dDsc_dt, dC_RIh_dt, dC_SIh_dt, dC_RIc_dt, dC_SIc_dt)))
  })
}

# --- 5. Solver Function ---
solve_model <- function(y_init, times, params) {
  out <- tryCatch({
    ode(y = y_init, times = times, func = infection_model_ode, parms = params,
        method = "radau", rtol = 1e-10, atol = 1e-12, hmax = 1)
  }, error = function(e) {
    message("radau failed, trying lsoda: ", e$message)
    ode(y = y_init, times = times, func = infection_model_ode, parms = params,
        method = "lsoda", rtol = 1e-8, atol = 1e-10, hmax = 1)
  })
  
  if(any(is.na(out))) {
    message("ODE solver produced NA values. Returning NULL.")
    return(NULL)
  }
  
  return(as.data.frame(out))
}

# --- 6. Function to Calculate Annual Outcomes (FINAL, DETAILED VERSION) ---
# --- THIS IS THE REQUIRED FUNCTION. 
calculate_annual_outcomes <- function(output_df, economic_start_year, economic_horizon_years) {
  if (is.null(output_df)) return(NULL)
  
  # Get the economic period subset
  econ_start_time <- economic_start_year * 365
  econ_end_time <- (economic_start_year + economic_horizon_years) * 365
  start_idx <- which(output_df$time >= econ_start_time)[1]
  end_idx <- which(output_df$time <= econ_end_time)
  end_idx <- end_idx[length(end_idx)]
  
  if (is.na(start_idx) || is.na(end_idx)) {
    return(NULL)
  }
  
  econ_data <- output_df[start_idx:end_idx, ]
  
  annual_outcomes <- list()
  for (i in 2:nrow(econ_data)) {
    year_idx <- i - 1
    year_name <- paste0("year_", year_idx)
    
    # ANNUAL DEATHS (Components and Totals)
    annual_res_deaths_h <- econ_data$Drh[i] - econ_data$Drh[i-1]
    annual_sus_deaths_h <- econ_data$Dsh[i] - econ_data$Dsh[i-1]
    annual_total_deaths_h <- annual_res_deaths_h + annual_sus_deaths_h
    
    annual_res_deaths_c <- econ_data$Drc[i] - econ_data$Drc[i-1]
    annual_sus_deaths_c <- econ_data$Dsc[i] - econ_data$Dsc[i-1]
    annual_total_deaths_c <- annual_res_deaths_c + annual_sus_deaths_c
    
    # ANNUAL NEW CASES (by difference from previous timepoint)
    annual_res_cases_h <- econ_data$C_RIh[i] - econ_data$C_RIh[i - 1]
    annual_sus_cases_h <- econ_data$C_SIh[i] - econ_data$C_SIh[i - 1]
    annual_total_cases_h <- annual_res_cases_h + annual_sus_cases_h
    
    annual_res_cases_c <- econ_data$C_RIc[i] - econ_data$C_RIc[i - 1]
    annual_sus_cases_c <- econ_data$C_SIc[i] - econ_data$C_SIc[i - 1]
    annual_total_cases_c <- annual_res_cases_c + annual_sus_cases_c
    
    # Return a list with the full breakdown
    annual_outcomes[[year_name]] <- list(
      annual_res_deaths_h = annual_res_deaths_h,
      annual_sus_deaths_h = annual_sus_deaths_h,
      annual_total_deaths_h = annual_total_deaths_h,
      annual_res_deaths_c = annual_res_deaths_c,
      annual_sus_deaths_c = annual_sus_deaths_c,
      annual_total_deaths_c = annual_total_deaths_c,
      
      annual_res_cases_h = annual_res_cases_h,
      annual_sus_cases_h = annual_sus_cases_h,
      annual_total_cases_h = annual_total_cases_h,
      annual_res_cases_c = annual_res_cases_c,
      annual_sus_cases_c = annual_sus_cases_c,
      annual_total_cases_c = annual_total_cases_c
    )
  }
  return(annual_outcomes)
}



# --- 7. Probabilistic Sensitivity Analysis (PSA) ---
set.seed(123)
n_psa_samples <- 1000

cat("\nRunning PSA with", n_psa_samples, "samples...\n")
cat("MCMC posterior has", nrow(posterior_samples_mcmc_selected), "samples available.\n")

start_time <- Sys.time()

# Run PSA simulations in parallel
psa_results_list <- foreach(i = 1:n_psa_samples, .packages = c("deSolve", "dplyr")) %dopar% {
  
  # 1. Randomly sample from MCMC posterior
  random_idx <- sample(nrow(posterior_samples_mcmc_selected), 1)
  mcmc_draw <- as.list(posterior_samples_mcmc_selected[random_idx, ])
  
  # 2. Build parameter sets for both scenarios
  params_baseline <- c(mcmc_draw, baseline_fixed_params)
  params_no_res <- params_baseline 
  params_no_res$delta_rch_rih <- 0
  params_no_res$delta_rcc_ric <- 0
  
  # 3. & 4. Run both scenarios
  out_baseline <- try(ode(y = y_init_baseline, times = times, func = infection_model_ode, parms = params_baseline, method = "lsoda"), silent = TRUE)
  out_no_res <- try(ode(y = y_init_no_resistance, times = times, func = infection_model_ode, parms = params_no_res, method = "lsoda"), silent = TRUE)
  
  # 5. Check for simulation failure
  if (inherits(out_baseline, "try-error") || inherits(out_no_res, "try-error")) {
    return(NULL)
  }
  
  # 6. Calculate annual outcomes using the NEW detailed function
  economic_start_year <- full_sim_years - economic_horizon_years
  baseline_outcomes <- calculate_annual_outcomes(as.data.frame(out_baseline), economic_start_year, economic_horizon_years)
  no_res_outcomes <- calculate_annual_outcomes(as.data.frame(out_no_res), economic_start_year, economic_horizon_years)
  
  if (is.null(baseline_outcomes) || is.null(no_res_outcomes)) return(NULL)
  
  # Prepare annual data frame. bind_rows will now create ALL the detailed columns.
  annual_baseline_df <- bind_rows(baseline_outcomes, .id = "year") %>% mutate(scenario = "baseline")
  annual_no_res_df <- bind_rows(no_res_outcomes, .id = "year") %>% mutate(scenario = "no_resistance")
  annual_results_df <- bind_rows(annual_baseline_df, annual_no_res_df)
  
  # Calculate TOTAL outcomes for the summary data frame
  total_baseline_deaths_h <- sum(sapply(baseline_outcomes, function(x) x$annual_total_deaths_h))
  total_baseline_deaths_c <- sum(sapply(baseline_outcomes, function(x) x$annual_total_deaths_c))
  
  total_no_res_deaths_h <- sum(sapply(no_res_outcomes, function(x) x$annual_total_deaths_h))
  total_no_res_deaths_c <- sum(sapply(no_res_outcomes, function(x) x$annual_total_deaths_c))
  
  # Create the summary data frame for this iteration
  summary_df <- data.frame(
    sample_id = i,
    incremental_deaths_h = total_baseline_deaths_h - total_no_res_deaths_h,
    incremental_deaths_c = total_baseline_deaths_c - total_no_res_deaths_c
  )
  
  return(list(summary = summary_df, annual = annual_results_df))
}



# Stop parallel processing
stopCluster(cl)

# --- Combine and Process All PSA Results ---
successful_results <- psa_results_list[!sapply(psa_results_list, is.null)]
psa_summary_results <- bind_rows(lapply(successful_results, function(x) x$summary))
psa_annual_results <- bind_rows(lapply(successful_results, function(x) x$annual), .id = "sample_id") %>%
  mutate(sample_id = as.integer(sample_id))

# --- Report Completion Status ---
n_successful <- nrow(psa_summary_results)
n_failed <- n_psa_samples - n_successful
total_time <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\PSA complete in %.1f minutes!\n", total_time))
cat(sprintf("Successfully completed: %d/%d iterations (%.1f%%)\n", 
            n_successful, n_psa_samples, (n_successful/n_psa_samples)*100))
if (n_failed > 0) {
  cat(sprintf("Failed iterations: %d\n", n_failed))
}


library(dplyr)

a_baseline_summary <- psa_annual_results %>%
  filter(scenario == "baseline") %>%
  group_by(year) %>%
  summarise(
    res_deaths_h_mean = mean(annual_res_deaths_h),
    res_deaths_h_median = median(annual_res_deaths_h),
    res_deaths_h_lower = quantile(annual_res_deaths_h, 0.025),
    res_deaths_h_upper = quantile(annual_res_deaths_h, 0.975),
    
    res_deaths_c_mean = mean(annual_res_deaths_c),
    res_deaths_c_median = median(annual_res_deaths_c),
    res_deaths_c_lower = quantile(annual_res_deaths_c, 0.025),
    res_deaths_c_upper = quantile(annual_res_deaths_c, 0.975),
    
    total_deaths_res_mean = mean(annual_res_deaths_h + annual_res_deaths_c),
    total_deaths_res_median = median(annual_res_deaths_h + annual_res_deaths_c),
    total_deaths_res_lower = quantile(annual_res_deaths_h + annual_res_deaths_c, 0.025),
    total_deaths_res_upper = quantile(annual_res_deaths_h + annual_res_deaths_c, 0.975),
    
    total_deaths_h_mean = mean(annual_sus_deaths_h),
    sus_deaths_h_median = median(annual_sus_deaths_h),
    sus_deaths_h_lower = quantile(annual_sus_deaths_h, 0.025),
    sus_deaths_h_upper = quantile(annual_sus_deaths_h, 0.975),
    
    sus_deaths_c_mean = mean(annual_sus_deaths_c),
    sus_deaths_c_median = median(annual_sus_deaths_c),
    sus_deaths_c_lower = quantile(annual_sus_deaths_c, 0.025),
    sus_deaths_c_upper = quantile(annual_sus_deaths_c, 0.975),
    
    total_deaths_h_mean = mean(annual_total_deaths_h),
    total_deaths_h_median = median(annual_total_deaths_h),
    total_deaths_h_lower = quantile(annual_total_deaths_h, 0.025),
    total_deaths_h_upper = quantile(annual_total_deaths_h, 0.975),
    
    total_deaths_c_mean = mean(annual_total_deaths_c),
    total_deaths_c_median = median(annual_total_deaths_c),
    total_deaths_c_lower = quantile(annual_total_deaths_c, 0.025),
    total_deaths_c_upper = quantile(annual_total_deaths_c, 0.975)
  )



# Filter to baseline scenario only
baseline_only <- psa_annual_results %>%
  filter(scenario == "baseline") %>%
  mutate(year = as.integer(gsub("year_", "", year)))

# Summarize cumulative resistant infections per year
annual_ri_summary <- baseline_only %>%
  group_by(year) %>%
  group_by(year) %>%
  summarise(
    RIh_mean = mean(annual_res_cases_h),
    RIh_median = median(annual_res_cases_h),
    RIh_lower = quantile(annual_res_cases_h, 0.025),
    RIh_upper = quantile(annual_res_cases_h, 0.975),
    
    RIc_mean = mean(annual_res_cases_c),
    RIc_median = median(annual_res_cases_c),
    RIc_lower = quantile(annual_res_cases_c, 0.025),
    RIc_upper = quantile(annual_res_cases_c, 0.975),
    
    RI_mean = mean(annual_res_cases_h + annual_res_cases_c),
    RI_median = median(annual_res_cases_h + annual_res_cases_c),
    RI_lower = quantile(annual_res_cases_h + annual_res_cases_c, 0.025),
    RI_upper = quantile(annual_res_cases_h + annual_res_cases_c, 0.975)
  )


