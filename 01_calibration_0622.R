## =============================================================================
## OHARP model calibration
##
## Clean appendix version of 01_calibraion_0620.R.
##
## Usage:
##   Rscript 01_calibration_0620_clean.R
##
## Requirements:
##   - model_calibration.stan must be in the same folder as this R script.
##   - Required R packages must already be installed.
##
## Optional environment overrides:
##   CALIBRATION_ITER, CALIBRATION_CHAINS, CALIBRATION_SEED,
##   CALIBRATION_ADAPT_DELTA, CALIBRATION_MAX_TREEDEPTH,
##   CALIBRATION_OUTPUT_DIR
## =============================================================================

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_idx <- grep("^--file=", cmd_args)
  if (length(file_idx) > 0) {
    script_path <- sub("^--file=", "", cmd_args[file_idx[1]])
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }

  file_flag_idx <- which(cmd_args == "--file")
  if (length(file_flag_idx) > 0 && file_flag_idx[1] < length(cmd_args)) {
    script_path <- cmd_args[file_flag_idx[1] + 1]
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }

  frame_files <- vapply(
    sys.frames(),
    function(frame) if (!is.null(frame$ofile)) frame$ofile else NA_character_,
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[length(frame_files)], mustWork = FALSE)))
  }

  getwd()
}

script_dir <- get_script_dir()
stan_file <- file.path(script_dir, "model_calibration.stan")
if (!file.exists(stan_file)) {
  stop("Cannot find model_calibration.stan in the script folder: ", script_dir)
}

output_dir <- Sys.getenv(
  "CALIBRATION_OUTPUT_DIR",
  unset = file.path(script_dir, "outputs_calibration_0620")
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_packages <- c("rstan", "bayesplot", "ggplot2", "dplyr", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Please install them before running this script."
  )
}

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)

rstan_options(auto_write = TRUE)
options(mc.cores = max(1L, parallel::detectCores()))

## -----------------------------------------------------------------------------
## 1. Initial conditions
## -----------------------------------------------------------------------------

N <- 5.18e6
Nh <- 37588 * 0.857
Nc <- N - Nh

RCh <- round(0.124 * Nh)
RIh <- round(0.119 * 0.314 * 0.377 * Nh)
SIh <- round(0.119 * 0.314 * 0.56 * Nh)
Sh <- Nh - RCh - RIh - SIh

RCc <- round(0.257 * Nc)
RIc <- round(0.001 * Nc)
SIc <- round(0.01 * Nc)
Sc <- Nc - RCc - RIc - SIc

y0 <- c(
  Sh = Sh,
  RCh = RCh,
  RIh = RIh,
  SIh = SIh,
  Sc = Sc,
  RCc = RCc,
  RIc = RIc,
  SIc = SIc,
  Drh = 0,
  Dsh = 0,
  Drc = 0,
  Dsc = 0,
  C_RIh = 0,
  C_SIh = 0,
  C_RIc = 0,
  C_SIc = 0
)

living_total <- Sh + RCh + RIh + SIh + Sc + RCc + RIc + SIc
if (!isTRUE(all.equal(living_total, N, tolerance = 1e-8))) {
  warning(
    "Initial living compartments sum to ",
    living_total,
    ", not total population ",
    N,
    "."
  )
}

## -----------------------------------------------------------------------------
## 2. Observed calibration targets, 2012-2021
## -----------------------------------------------------------------------------

incidence_per_10k <- c(38, 35, 34.5, 33.5, 37, 27.5, 28.5, 28.2, 27, 24.5)
yearly_data <- data.frame(
  year = 2012:2021,
  RIh_obs = round(incidence_per_10k)
)

RIh_obs <- yearly_data$RIh_obs
n_years <- nrow(yearly_data)
n_days <- 365 * n_years
t0 <- 0
ts <- seq(1, n_days, by = 1)

## -----------------------------------------------------------------------------
## 3. Fixed model parameters
## -----------------------------------------------------------------------------

mu_sih <- 0.00311667
mu_rih <- 0.00404500
delta_rch_rih <- 0.03395087

gamma_sih_sh <- 0.142857
gamma_rch_sh <- 0.00612
gamma_rih_rch <- 0.1111

beta_sc_sic <- 0.00000792

gamma_sic_sc <- 0.25
gamma_rcc_sc <- 0.009
gamma_ric_rcc <- 0.17

mu_sic <- 0.003226
mu_ric <- 0.003629

alpha_adm <- 0.000316
alpha_dis <- 0.235

mu_b <- 0.00002518
mu_m <- 0.000002
mu_d <- 0.00001353

delta_rch_sih <- 0.00084658
delta_rcc_sic <- 0.00084658

data_sir <- list(
  n_years = n_years,
  n_days = n_days,
  y0 = y0,
  t0 = t0,
  ts = ts,
  RIh_obs = RIh_obs,

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

saveRDS(data_sir, file.path(output_dir, "calibration_data_sir.rds"))
write.csv(yearly_data, file.path(output_dir, "calibration_observed_incidence.csv"), row.names = FALSE)

## -----------------------------------------------------------------------------
## 4. Run Stan calibration
## -----------------------------------------------------------------------------

calibration_iter <- as.integer(Sys.getenv("CALIBRATION_ITER", unset = "4000"))
calibration_chains <- as.integer(Sys.getenv("CALIBRATION_CHAINS", unset = "4"))
calibration_seed <- as.integer(Sys.getenv("CALIBRATION_SEED", unset = "123"))
adapt_delta <- as.numeric(Sys.getenv("CALIBRATION_ADAPT_DELTA", unset = "0.99"))
max_treedepth <- as.integer(Sys.getenv("CALIBRATION_MAX_TREEDEPTH", unset = "15"))

if (is.na(calibration_iter) || calibration_iter <= 0) {
  stop("CALIBRATION_ITER must be a positive integer.")
}
if (is.na(calibration_chains) || calibration_chains <= 0) {
  stop("CALIBRATION_CHAINS must be a positive integer.")
}
if (is.na(calibration_seed)) {
  stop("CALIBRATION_SEED must be an integer.")
}
if (is.na(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1) {
  stop("CALIBRATION_ADAPT_DELTA must be between 0 and 1.")
}
if (is.na(max_treedepth) || max_treedepth <= 0) {
  stop("CALIBRATION_MAX_TREEDEPTH must be a positive integer.")
}

cat("Stan model:", stan_file, "\n")
cat("Output directory:", output_dir, "\n")
cat(
  "Running calibration with iter =",
  calibration_iter,
  ", chains =",
  calibration_chains,
  ", seed =",
  calibration_seed,
  "\n"
)

fit6 <- rstan::stan(
  file = stan_file,
  data = data_sir,
  iter = calibration_iter,
  chains = calibration_chains,
  seed = calibration_seed,
  control = list(
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
)

saveRDS(fit6, file.path(output_dir, "fit6.rds"))
save(fit6, file = file.path(output_dir, "fit6.RData"))

## -----------------------------------------------------------------------------
## 5. Diagnostics and posterior outputs
## -----------------------------------------------------------------------------

pars_calibrated <- c(
  "beta_hh",
  "beta_eh",
  "beta_hc",
  "beta_sc_ac",
  "beta_sc_ec",
  "beta_sh_sih",
  "delta_rcc_ric"
)

capture.output(
  print(fit6, digits = 9, pars = pars_calibrated),
  file = file.path(output_dir, "calibrated_parameter_summary.txt")
)

posterior_draws <- as.data.frame(fit6, pars = pars_calibrated)
saveRDS(posterior_draws, file.path(output_dir, "calibrated_parameter_draws.rds"))
write.csv(posterior_draws, file.path(output_dir, "calibrated_parameter_draws.csv"), row.names = FALSE)

posterior_summary <- posterior_draws %>%
  summarise(across(
    everything(),
    list(
      mean = ~ mean(.x, na.rm = TRUE),
      sd = ~ sd(.x, na.rm = TRUE),
      q025 = ~ quantile(.x, 0.025, na.rm = TRUE),
      median = ~ median(.x, na.rm = TRUE),
      q975 = ~ quantile(.x, 0.975, na.rm = TRUE)
    )
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("parameter", ".value"),
    names_pattern = "(.+)_(mean|sd|q025|median|q975)$"
  )
write.csv(posterior_summary, file.path(output_dir, "calibrated_parameter_summary.csv"), row.names = FALSE)

minmax <- t(vapply(
  posterior_draws,
  function(x) c(min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE)),
  c(min = 0, max = 0)
))
write.csv(
  round(as.data.frame(minmax), 9),
  file.path(output_dir, "calibrated_parameter_minmax.csv")
)

p_trace <- rstan::traceplot(fit6, pars = pars_calibrated)
ggplot2::ggsave(
  filename = file.path(output_dir, "calibrated_parameter_traceplot.png"),
  plot = p_trace,
  width = 12,
  height = 8,
  dpi = 300
)

p_density <- rstan::stan_dens(fit6, pars = pars_calibrated, separate_chains = TRUE)
ggplot2::ggsave(
  filename = file.path(output_dir, "calibrated_parameter_density.png"),
  plot = p_density,
  width = 12,
  height = 8,
  dpi = 300
)

## -----------------------------------------------------------------------------
## 6. Posterior predictive check
## -----------------------------------------------------------------------------

incidence_post <- as.matrix(fit6, pars = "pred_yearly_total_RIh_incidence")
obs_incidence <- yearly_data$RIh_obs

write.csv(
  as.data.frame(incidence_post),
  file.path(output_dir, "posterior_pred_yearly_total_RIh_incidence.csv"),
  row.names = FALSE
)

bayesplot::color_scheme_set("brightblue")
p_ppc <- bayesplot::ppc_ribbon(
  y = obs_incidence,
  yrep = incidence_post,
  x = yearly_data$year,
  prob = 0.5,
  prob_outer = 0.95
) +
  geom_point(
    data = yearly_data,
    aes(x = year, y = RIh_obs),
    color = "navy",
    size = 3,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Year",
    y = "Incidence per 10,000 inpatient days",
    title = "Posterior Predictive Intervals vs Observed Incidence",
    subtitle = "Showing 50% and 95% posterior predictive intervals"
  ) +
  theme_minimal()

ggplot2::ggsave(
  filename = file.path(output_dir, "ppc_ribbon_plot.png"),
  plot = p_ppc,
  width = 10,
  height = 6,
  dpi = 300
)

cat("\nCalibration complete. Outputs saved to:\n")
cat(output_dir, "\n")
