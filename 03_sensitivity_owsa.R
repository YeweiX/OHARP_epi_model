# ==============================================================================
# Analysis: One-Way Sensitivity Analysis (OWSA)
# Description: Assesses impact of parameters on 10-Year Cumulative Burden.
# Method: 1. Rank parameters (Tornado) -> 2. Detailed Sweep -> 3. Spider Plot
# ==============================================================================

# --- 1. SETUP & INITIALIZATION ------------------------------------------------
rm(list = ls())

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

# 1.1. Define Baseline Parameters (Calibrated Values)
baseline_params <- list(
  # Transmission
  beta_hh = 0.007893581, beta_hc = 0.003655865,
  beta_eh = 0.006259864, beta_sc_ec = 0.000260089, beta_sc_ac = 0.000683476,
  beta_sh_sih = 0.239527636, beta_sc_sic = 0.00000792,
  
  # Progression
  delta_rch_rih = 0.03395087, delta_rch_sih = 0.00084658,
  delta_rcc_ric = 0.000085556, delta_rcc_sic = 0.00084658,
  
  # Recovery & Mortality
  mu_sih = 0.00311667, mu_rih = 0.00404500, 
  gamma_sih_sh = 0.14285714, gamma_rch_sh = 0.00612032, gamma_rih_rch = 0.11111111,
  gamma_sic_sc = 0.25, gamma_rcc_sc = 0.009, gamma_ric_rcc = 0.17,
  mu_sic = 0.003226, mu_ric = 0.003629,
  
  # Demographics & Movement
  alpha_adm = 0.000316, alpha_dis = 0.235,
  mu_b = 0.00002518, mu_m = 0.000002, mu_d = 0.00001353
)

# 1.2. Define Initial Conditions (Matches Scenario Analysis)
N_total <- 5.18e6
N_hosp  <- 37588 * 0.857
N_comm  <- N_total - N_hosp

# Hospital Compartments
RC_h <- round(0.124 * N_hosp)
RI_h <- round(0.119 * 0.314 * 0.377 * N_hosp)
SI_h <- round(0.119 * 0.314 * 0.56 * N_hosp)
S_h  <- N_hosp - RC_h - RI_h - SI_h

# Community Compartments
RC_c <- round(0.257 * N_comm)
RI_c <- round(0.001 * N_comm)
SI_c <- round(0.01 * N_comm)
S_c  <- N_comm - RC_c - RI_c - SI_c

y_init_baseline <- c(
  Sh = S_h, RCh = RC_h, RIh = RI_h, SIh = SI_h,
  Sc = S_c, RCc = RC_c, RIc = RI_c, SIc = SI_c,
  Drh = 0, Dsh = 0, Drc = 0, Dsc = 0,        # Cumulative Deaths
  C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0 # Cumulative Cases
)

# --- 2. MODEL DEFINITION (Self-Contained) -------------------------------------

# Helper: Constrain flow to avoid negative populations
constrain_flow <- function(rate, pool) { min(rate, pool) }

# ODE System
infection_model_ode <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Populations
    Nh_curr <- max(1, Sh + RCh + RIh + SIh)
    Nc_curr <- max(1, Sc + RCc + RIc + SIc)
    N_living <- Nh_curr + Nc_curr
    
    # --- Flows: Hospital ---
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
    
    # --- Flows: Community ---
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

# --- 3. SIMULATION WRAPPER ----------------------------------------------------

# Settings
full_sim_years <- 20
econ_years <- 10
start_day_econ <- (full_sim_years - econ_years) * 365
times_vec <- seq(0, full_sim_years * 365, by = 365) # Yearly steps

run_simulation_horizon <- function(parms) {
  # Run ODE
  out <- try(ode(y = y_init_baseline, times = times_vec, 
                 func = infection_model_ode, parms = parms, method = "lsoda"), 
             silent = TRUE)
  
  if (inherits(out, "try-error")) return(rep(NA, 4))
  
  df <- as.data.frame(out)
  
  # Extract start and end of economic horizon
  row_start <- df[df$time >= start_day_econ, ][1, ]
  row_end   <- df[nrow(df), ]
  
  if (any(is.na(row_start)) || any(is.na(row_end))) return(rep(NA, 4))
  
  # Calculate 10-Year Cumulative Totals (End - Start)
  c(
    Drh = row_end$Drh - row_start$Drh,
    Drc = row_end$Drc - row_start$Drc,
    Cases_h = row_end$C_RIh - row_start$C_RIh,
    Cases_c = row_end$C_RIc - row_start$C_RIc
  )
}

# --- 4. SENSITIVITY ANALYSIS --------------------------------------------------

# 4.1. Pre-Scan: Identify Top Influential Parameters
message("Step 1: Identifying top influential parameters...")

# Structural parameters to exclude from sensitivity
exclude_pars <- c("alpha_adm", "alpha_dis", "mu_b", "mu_m", "mu_d")
scan_results <- list()

for (p_name in names(baseline_params)) {
  if (p_name %in% exclude_pars) next
  val <- baseline_params[[p_name]]
  if (val == 0) next
  
  # Check +/- 50% range
  res_low  <- run_simulation_horizon(replace(baseline_params, p_name, val * 0.5))
  res_high <- run_simulation_horizon(replace(baseline_params, p_name, val * 1.5))
  
  if (!any(is.na(res_low)) && !any(is.na(res_high))) {
    scan_results[[p_name]] <- data.frame(
      parameter = p_name,
      range_Drh = abs(res_high["Drh"] - res_low["Drh"]),
      range_Cases_c = abs(res_high["Cases_c"] - res_low["Cases_c"])
    )
  }
}

df_scan <- bind_rows(scan_results)

# Select Top 5 for Hospital Deaths and Community Cases
top5_Drh <- df_scan %>% arrange(desc(range_Drh)) %>% slice(1:5) %>% pull(parameter)
top5_Cases_c <- df_scan %>% arrange(desc(range_Cases_c)) %>% slice(1:5) %>% pull(parameter)
key_params <- unique(c(top5_Drh, top5_Cases_c))

message("Key parameters identified: ", paste(key_params, collapse=", "))

# 4.2. Detailed Sweep for Top Parameters
message("Step 2: Running detailed parameter sweep...")
pct_changes <- seq(-0.5, 0.5, by = 0.05) # +/- 50% sweep
baseline_res <- run_simulation_horizon(baseline_params)
sweep_data <- list()

for (p_name in key_params) {
  for (pct in pct_changes) {
    curr_pars <- baseline_params
    curr_pars[[p_name]] <- curr_pars[[p_name]] * (1 + pct)
    
    res <- run_simulation_horizon(curr_pars)
    
    if (!any(is.na(res))) {
      sweep_data[[length(sweep_data) + 1]] <- data.frame(
        parameter = p_name,
        pct_change = pct,
        # Calculate Absolute Change relative to Baseline
        delta_Drh = res["Drh"] - baseline_res["Drh"],
        delta_Cases_c = res["Cases_c"] - baseline_res["Cases_c"]
      )
    }
  }
}

df_sweep <- bind_rows(sweep_data)

# --- 5. VISUALIZATION (Spider Plot) -------------------------------------------

# Label Dictionary for Clean Plots
param_labels <- c(
  "beta_hh" = "Hosp. Trans (Human)", "beta_hc" = "Comm. Trans (Human)",
  "beta_sh_sih" = "Hosp. Trans (Sensitive)", "delta_rch_rih" = "Hosp. Progression",
  "mu_rih" = "Hosp. Mortality (Res)", "gamma_rih_rch" = "Hosp. Recovery",
  "beta_sc_ec" = "Env-to-Comm Trans", "beta_sc_ac" = "Animal-to-Comm Trans"
)

# Prepare Data for Plotting
df_plot <- df_sweep %>%
  mutate(label = recode(parameter, !!!param_labels, .default = parameter))

# Plot A: Sensitivity of Hospital Deaths
p1 <- df_plot %>%
  filter(parameter %in% top5_Drh) %>%
  ggplot(aes(x = pct_change, y = delta_Drh, color = label, group = parameter)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  scale_x_continuous(labels = percent_format()) +
  labs(title = "Impact on Hospital Resistant Deaths",
       subtitle = "Change in 10-year cumulative deaths",
       y = "Change in Deaths (n)", x = "Parameter Variation (%)", color = "Parameter") +
  theme_bw()

# Plot B: Sensitivity of Community Cases
p2 <- df_plot %>%
  filter(parameter %in% top5_Cases_c) %>%
  ggplot(aes(x = pct_change, y = delta_Cases_c, color = label, group = parameter)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  scale_x_continuous(labels = percent_format()) +
  labs(title = "Impact on Community Resistant Infections",
       subtitle = "Change in 10-year cumulative infections",
       y = "Change in Cases (n)", x = "Parameter Variation (%)", color = "Parameter") +
  theme_bw()

# Combine
final_plot <- p1 / p2 + plot_layout(guides = "collect")
print(final_plot)

ggsave("Figure_Sensitivity_Spider.png", final_plot, width = 8, height = 10)
