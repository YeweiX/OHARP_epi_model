# ==============================================================================
# Analysis: Cumulative Burden of Resistance (Counterfactual Analysis)
# Context: Hospital (h) vs. Community (c) Transmission Model
# Methodology: Probabilistic Sensitivity Analysis (PSA) with 10-Year Horizon
# Output: Cumulative attributable mortality/morbidity with 95% Credible Intervals
# ==============================================================================

# --- 1. Setup and Environment -------------------------------------------------
rm(list = ls())

# Load required packages
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(doParallel)
library(foreach)

# Setup Parallel Processing
num_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Parallel processing initialized on", num_cores, "cores.\n")

# Load Fitted Model Data
# NOTE: Ensure 'fit_esbl_results.RData' (or fit6.RData) is in your working directory
if(file.exists("fit_esbl_results.RData")) {
  load("fit_esbl_results.RData")
} else {
  warning("Model file not found. Please load the posterior samples manually.")
}

# --- 2. Parameter Definitions -------------------------------------------------

# Fixed Parameters (Biological & Demographic Baseline)
# _h = Hospital, _c = Community
fixed_params <- list(
  # Transmission & Progression
  beta_hh = 0.00789, beta_sh_sih = 0.23953, 
  delta_rch_rih = 0.03395, delta_rch_sih = 0.00085,
  beta_sc_sic = 7.92e-6, delta_rcc_ric = 8.56e-5, delta_rcc_sic = 0.00085,
  beta_eh = 0.00626, beta_hc = 0.00366, beta_sc_ac = 0.00068, beta_sc_ec = 0.00026,
  
  # Mortality & Recovery
  mu_sih = 0.00312, mu_rih = 0.00405, 
  gamma_sih_sh = 0.14286, gamma_rch_sh = 0.00612, gamma_rih_rch = 0.11111,
  gamma_sic_sc = 0.25, gamma_rcc_sc = 0.009, gamma_ric_rcc = 0.17,
  mu_sic = 0.00323, mu_ric = 0.00363,
  
  # Demographics
  alpha_adm = 0.00032, alpha_dis = 0.235,
  mu_b = 2.52e-5, mu_m = 2e-6, mu_d = 1.35e-5
)

# Extract MCMC Samples for PSA
mcmc_pars <- c("beta_hh", "beta_eh", "beta_hc", "beta_sc_ac", "beta_sc_ec", "beta_sh_sih", "delta_rcc_ric")
# Assuming 'fit6' is the loaded stan object
if(exists("fit6")) {
  posterior_samples <- as.data.frame(rstan::extract(fit6, pars = mcmc_pars))
} else {
  stop("Stan fit object 'fit6' not found.")
}

# --- 3. Model Logic -----------------------------------------------------------

# ODE Function
ode_model <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Populations
    Nh <- max(1e-6, Sh + RCh + RIh + SIh)
    Nc <- max(1e-6, Sc + RCc + RIc + SIc)
    
    # Hospital Flows
    # ... (Full differential equations abbreviated for brevity, same as previous scripts)
    # Key concept: We track cumulative deaths (Drh, Drc) to calculate incidence later
    
    # Re-using the standard equations defined in previous steps...
    # (Ensure the full ODE function from the previous response is included here in the actual run)
    list(c(dSh, dRCh, dRIh, dSIh, dSc, dRCc, dRIc, dSIc,
           dDrh, dDsh, dDrc, dDsc, dC_RIh, dC_SIh, dC_RIc, dC_SIc))
  })
}
# Note: For actual execution, insert the full `ode_model` function body here.

# Helper: Calculate Annual Outcomes from Cumulative Time Series
calc_annual_stats <- function(df, start_yr, horizon_yrs) {
  # Filter time window
  t_start <- start_yr * 365
  t_end   <- (start_yr + horizon_yrs) * 365
  df_sub  <- df[df$time >= t_start & df$time <= t_end, ]
  
  # Calculate differences (Year t - Year t-1)
  stats_list <- list()
  for(i in 2:nrow(df_sub)) {
    year_num <- i - 1
    
    # Deaths
    deaths_h_res <- df_sub$Drh[i] - df_sub$Drh[i-1]
    deaths_h_sus <- df_sub$Dsh[i] - df_sub$Dsh[i-1]
    deaths_c_res <- df_sub$Drc[i] - df_sub$Drc[i-1]
    deaths_c_sus <- df_sub$Dsc[i] - df_sub$Dsc[i-1]
    
    # Cases
    cases_h_res <- df_sub$C_RIh[i] - df_sub$C_RIh[i-1]
    cases_h_sus <- df_sub$C_SIh[i] - df_sub$C_SIh[i-1]
    cases_c_res <- df_sub$C_RIc[i] - df_sub$C_RIc[i-1]
    cases_c_sus <- df_sub$C_SIc[i] - df_sub$C_SIc[i-1]
    
    stats_list[[i-1]] <- data.frame(
      year = year_num,
      deaths_h_total = deaths_h_res + deaths_h_sus,
      deaths_c_total = deaths_c_res + deaths_c_sus,
      deaths_h_res   = deaths_h_res,
      deaths_c_res   = deaths_c_res,
      cases_h_total  = cases_h_res + cases_h_sus,
      cases_c_total  = cases_c_res + cases_c_sus
    )
  }
  bind_rows(stats_list)
}

# --- 4. PSA Simulation Loop ---------------------------------------------------
n_psa <- 1000
years_sim <- 20
years_econ <- 10
times <- seq(0, years_sim * 365, by = 365)

# Initial Conditions (Baseline)
N_total <- 5.18e6; N_h <- 37588 * 0.857; N_c <- N_total - N_h
y_base <- c(Sh=N_h*0.4, RCh=N_h*0.12, RIh=N_h*0.01, SIh=N_h*0.02, # Approx
            Sc=N_c*0.7, RCc=N_c*0.25, RIc=N_c*0.001, SIc=N_c*0.01,
            Drh=0, Dsh=0, Drc=0, Dsc=0, C_RIh=0, C_SIh=0, C_RIc=0, C_SIc=0)

# Initial Conditions (No Resistance: R -> S)
y_nores <- y_base
y_nores["Sh"] <- y_base["Sh"] + y_base["RCh"] + y_base["RIh"]
y_nores["Sc"] <- y_base["Sc"] + y_base["RCc"] + y_base["RIc"]
y_nores[c("RCh","RIh","RCc","RIc")] <- 0

message("Starting PSA...")

psa_output <- foreach(i = 1:n_psa, .packages = c("deSolve", "dplyr")) %dopar% {
  
  # 1. Parameter Sampling
  sample_idx <- sample(nrow(posterior_samples), 1)
  pars_curr  <- c(as.list(posterior_samples[sample_idx,]), fixed_params)
  
  # 2. Define Counterfactual Parameters (Transmission to R = 0)
  pars_nores <- pars_curr
  pars_nores$delta_rch_rih <- 0
  pars_nores$delta_rcc_ric <- 0
  
  # 3. Run Models
  out_base  <- as.data.frame(ode(y=y_base, times=times, func=ode_model, parms=pars_curr))
  out_nores <- as.data.frame(ode(y=y_nores, times=times, func=ode_model, parms=pars_nores))
  
  # 4. Calculate Annual Stats
  start_yr <- years_sim - years_econ
  stats_base  <- calc_annual_stats(out_base, start_yr, years_econ)
  stats_nores <- calc_annual_stats(out_nores, start_yr, years_econ)
  
  # 5. Calculate Cumulative & Incremental Totals (The Key Step)
  # Sum over the 10-year horizon
  cum_base <- colSums(stats_base[, -1]) # Exclude year column
  cum_nores <- colSums(stats_nores[, -1])
  
  # Incremental Burden = Baseline - No_Resistance
  # (Positive number = Deaths caused by resistance)
  burden_deaths_h <- cum_base["deaths_h_total"] - cum_nores["deaths_h_total"]
  burden_deaths_c <- cum_base["deaths_c_total"] - cum_nores["deaths_c_total"]
  
  # Return summarized row
  data.frame(
    iter = i,
    burden_deaths_h = burden_deaths_h,
    burden_deaths_c = burden_deaths_c,
    burden_deaths_total = burden_deaths_h + burden_deaths_c
  )
}

results_df <- bind_rows(psa_output)
stopCluster(cl)

# --- 5. Analysis & Visualization ----------------------------------------------

# Summary Statistics
summary_table <- results_df %>%
  summarise(
    mean_total = mean(burden_deaths_total),
    ci_low_total = quantile(burden_deaths_total, 0.025),
    ci_high_total = quantile(burden_deaths_total, 0.975),
    
    mean_h = mean(burden_deaths_h),
    mean_c = mean(burden_deaths_c)
  )

print(summary_table)

# Visualization: Boxplot of Cumulative Burden
plot_data <- results_df %>%
  pivot_longer(cols = starts_with("burden"), names_to = "Metric", values_to = "Deaths") %>%
  mutate(Metric = case_when(
    Metric == "burden_deaths_total" ~ "Overall",
    Metric == "burden_deaths_h" ~ "Hospital",
    Metric == "burden_deaths_c" ~ "Community"
  )) %>%
  mutate(Metric = factor(Metric, levels = c("Overall", "Hospital", "Community")))

p <- ggplot(plot_data, aes(x = Metric, y = Deaths, fill = Metric)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.2) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("#333333", "#0072B2", "#D55E00")) +
  labs(
    title = "Cumulative 10-Year Burden of Antibiotic Resistance",
    subtitle = "Attributable deaths (Baseline - No Resistance Counterfactual)",
    y = "Total Attributable Deaths",
    x = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 14))

print(p)
ggsave("Figure_Cumulative_Burden_Boxplot.png", p, width = 8, height = 6)
