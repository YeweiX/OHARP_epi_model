
// =============================================================================
// Stan Model: ESBL Transmission Dynamics (Open Cohort)
// Description: Two-setting (Hospital/Community) compartmental model 
//              calibrated to incidence data.
// =============================================================================

functions {
  
  // SIR-type ODE System
  real[] sir(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    
    // --- 1. Unpack State Variables ---
    // S: Susceptible, RC: Resistant Colonized, RI: Resistant Infected, SI: Sensitive Infected
    // Suffix: _h (Hospital), _c (Community)
    real Sh  = y[1];
    real RCh = y[2];
    real RIh = y[3];
    real SIh = y[4];
    
    real Sc  = y[5];
    real RCc = y[6];
    real RIc = y[7];
    real SIc = y[8];
    
    // Mortality tracking (not used in feedback loop, but tracked)
    // Dr: Death Resistant, Ds: Death Sensitive
    real Drh = y[9];
    real Dsh = y[10];
    real Drc = y[11];
    real Dsc = y[12];
    
    // Cumulative Incidence counters
    real C_RIh_total = y[13];
    real C_SIh_total = y[14];
    real C_RIc_total = y[15];
    real C_SIc_total = y[16];

    // --- 2. Unpack Estimated Parameters (theta) ---
    real beta_eh       = theta[1]; // Environmental-to-Hospital transmission
    real beta_sc_ec    = theta[2]; // Environmental-to-Community transmission
    real beta_sc_ac    = theta[3]; // Animal-to-Community transmission
    real beta_hh       = theta[4]; // Human-to-Human (Hospital)
    real beta_hc       = theta[5]; // Human-to-Human (Community)
    real beta_sh_sih   = theta[6]; // Sensitive Infection transmission (Hospital)
    real delta_rcc_ric = theta[7]; // Progression RC -> RI (Community)

    // --- 3. Unpack Fixed Parameters (x_r) ---
    real mu_sih        = x_r[1];
    real mu_rih        = x_r[2];
    real gamma_sih_sh  = x_r[3];
    real gamma_rch_sh  = x_r[4];
    real gamma_rih_rch = x_r[5];
    real beta_sc_sic   = x_r[6];
    real mu_sic        = x_r[7];
    real mu_ric        = x_r[8];
    real gamma_sic_sc  = x_r[9];
    real gamma_rcc_sc  = x_r[10];
    real gamma_ric_rcc = x_r[11];
    real delta_rch_sih = x_r[12];
    real delta_rch_rih = x_r[13];
    real alpha_adm     = x_r[14]; // Admission rate
    real alpha_dis     = x_r[15]; // Discharge rate
    real mu_b          = x_r[16]; // Birth rate
    real mu_m          = x_r[17]; // Migration rate
    real mu_d          = x_r[18]; // Natural death rate
    real delta_rcc_sic = x_r[19];

    // --- 4. Population Denominators ---
    real Nh_current     = fmax(1e-6, Sh + RCh + RIh + SIh);
    real Nc_current     = fmax(1e-6, Sc + RCc + RIc + SIc);
    real N_total_living = Nh_current + Nc_current;

    // --- 5. Differential Equations ---

    // Hospital: Susceptible (Sh)
    real dSh_dt = 
      - beta_eh * Sh                                    // Env transmission
      - beta_hh * Sh * RCh / Nh_current                 // Human transmission (RC)
      - beta_sh_sih * Sh * SIh / Nh_current             // Human transmission (SI)
      - alpha_dis * Sh                                  // Discharge
      - mu_d * Sh                                       // Death
      + gamma_sih_sh * SIh                              // Recovery
      + gamma_rch_sh * RCh                              // Clearance
      + alpha_adm * Sc;                                 // Admission

    // Hospital: Resistant Colonized (RCh)
    real dRCh_dt = 
      - delta_rch_rih * RCh                             // Progression to Infection
      - gamma_rch_sh * RCh                              // Clearance
      - alpha_dis * RCh                                 // Discharge
      - delta_rch_sih * RCh                             // Progression to Sensitive
      - mu_d * RCh                                      // Death
      + beta_eh * Sh                                    // Env Acquisition
      + beta_hh * Sh * RCh / Nh_current                 // Human Acquisition
      + gamma_rih_rch * RIh                             // Recovery from Infection
      + alpha_adm * RCc;                                // Admission

    // Hospital: Resistant Infected (RIh)
    real dRIh_dt = 
      - gamma_rih_rch * RIh                             // Recovery
      - mu_rih * RIh                                    // Mortality (Infection)
      + delta_rch_rih * RCh                             // Progression
      + alpha_adm * RIc;                                // Admission

    // Hospital: Sensitive Infected (SIh)
    real dSIh_dt = 
      - gamma_sih_sh * SIh                              // Recovery
      - mu_sih * SIh                                    // Mortality (Infection)
      + beta_sh_sih * Sh * SIh / Nh_current             // Transmission
      + delta_rch_sih * RCh                             // Transition
      + alpha_adm * SIc;                                // Admission

    // Community: Susceptible (Sc)
    real dSc_dt = 
      - beta_sc_ec * Sc 
      - beta_sc_ac * Sc * 0.33 
      - beta_hc * Sc * RCc / Nc_current 
      - beta_sc_sic * Sc * SIc / Nc_current 
      - alpha_adm * Sc 
      - mu_d * Sc 
      + gamma_sic_sc * SIc 
      + gamma_rcc_sc * RCc 
      + alpha_dis * Sh 
      + mu_b * N_total_living 
      + mu_m * N_total_living;

    // Community: Resistant Colonized (RCc)
    real dRCc_dt = 
      - delta_rcc_sic * RCc 
      - delta_rcc_ric * RCc 
      - gamma_rcc_sc * RCc 
      - alpha_adm * RCc 
      - mu_d * RCc 
      + beta_sc_ec * Sc 
      + beta_sc_ac * Sc * 0.33 
      + beta_hc * Sc * RCc / Nc_current 
      + gamma_ric_rcc * RIc 
      + alpha_dis * RCh;

    // Community: Resistant Infected (RIc)
    real dRIc_dt = 
      - gamma_ric_rcc * RIc 
      - mu_ric * RIc 
      - alpha_adm * RIc 
      + delta_rcc_ric * RCc;

    // Community: Sensitive Infected (SIc)
    real dSIc_dt = 
      - gamma_sic_sc * SIc 
      - mu_sic * SIc 
      - alpha_adm * SIc 
      + beta_sc_sic * Sc * SIc / Nc_current
      + delta_rcc_sic * RCc;

    // Mortality Accumulators
    real dDrh_dt = mu_rih * RIh;
    real dDsh_dt = mu_sih * SIh;
    real dDrc_dt = mu_ric * RIc;
    real dDsc_dt = mu_sic * SIc;

    // Cumulative Incidence Accumulators
    real dC_RIh_total_dt = delta_rch_rih * RCh + alpha_adm * RIc;
    real dC_SIh_total_dt = beta_sh_sih * Sh * SIh / Nh_current + delta_rch_sih * RCh + alpha_adm * SIc;
    real dC_RIc_total_dt = delta_rcc_ric * RCc;
    real dC_SIc_total_dt = beta_sc_sic * Sc * SIc / Nc_current + delta_rcc_sic * RCc;

    return {
      dSh_dt, dRCh_dt, dRIh_dt, dSIh_dt,
      dSc_dt, dRCc_dt, dRIc_dt, dSIc_dt,
      dDrh_dt, dDsh_dt, dDrc_dt, dDsc_dt,
      dC_RIh_total_dt, dC_SIh_total_dt, dC_RIc_total_dt, dC_SIc_total_dt
    };
  }
}

data {
  // Dimensions
  int<lower=1> n_years;
  int<lower=1> n_days;
  
  // Data Vectors
  real ts[n_days];         // Time steps
  real y0[16];             // Initial conditions
  real t0;                 // Start time
  int<lower=0> RIh_obs[n_years]; // Observed Incidence (Hospital Resistant)

  // Fixed Parameters (Clinical & Demographic)
  real<lower=0> mu_sih_fixed;
  real<lower=0> mu_rih_fixed;
  real<lower=0> gamma_sih_sh_fixed;
  real<lower=0> gamma_rch_sh_fixed;
  real<lower=0> gamma_rih_rch_fixed;
  real<lower=0> beta_sc_sic_fixed;
  real<lower=0> mu_sic_fixed;
  real<lower=0> mu_ric_fixed;
  real<lower=0> gamma_sic_sc_fixed;
  real<lower=0> gamma_rcc_sc_fixed;
  real<lower=0> gamma_ric_rcc_fixed;
  real<lower=0> delta_rch_sih_fixed;
  real<lower=0> delta_rch_rih_fixed;
  real<lower=0> alpha_adm_fixed;
  real<lower=0> alpha_dis_fixed;
  real<lower=0> mu_b_fixed;
  real<lower=0> mu_m_fixed;
  real<lower=0> mu_d_fixed;
  real<lower=0> delta_rcc_sic_fixed;
}

transformed data {
  // Map fixed parameters to x_r array for the ODE solver
  real x_r[19];
  x_r[1]  = mu_sih_fixed;
  x_r[2]  = mu_rih_fixed;
  x_r[3]  = gamma_sih_sh_fixed;
  x_r[4]  = gamma_rch_sh_fixed;
  x_r[5]  = gamma_rih_rch_fixed;
  x_r[6]  = beta_sc_sic_fixed;
  x_r[7]  = mu_sic_fixed;
  x_r[8]  = mu_ric_fixed;
  x_r[9]  = gamma_sic_sc_fixed;
  x_r[10] = gamma_rcc_sc_fixed;
  x_r[11] = gamma_ric_rcc_fixed;
  x_r[12] = delta_rch_sih_fixed;
  x_r[13] = delta_rch_rih_fixed;
  x_r[14] = alpha_adm_fixed;
  x_r[15] = alpha_dis_fixed;
  x_r[16] = mu_b_fixed;
  x_r[17] = mu_m_fixed;
  x_r[18] = mu_d_fixed;
  x_r[19] = delta_rcc_sic_fixed;

  int x_i[0]; // Empty integer array (not used but required by signature)
}

parameters {
  // Transmission Probabilities
  real<lower=0, upper=1>      beta_eh;       // Env -> Hosp
  real<lower=0, upper=1>      beta_sc_ec;    // Env -> Comm
  real<lower=0, upper=1>      beta_sc_ac;    // Animal -> Comm
  real<lower=0.0021, upper=0.0139> beta_hh;  // Human -> Hosp (Constrained)
  real<lower=0.0016, upper=0.011>  beta_hc;  // Human -> Comm (Constrained)
  real<lower=0, upper=1>      beta_sh_sih;   // Sensitive Trans. Hosp
  real<lower=0, upper=1>      delta_rcc_ric; // Progression Rate
  
  // Dispersion Parameter for NegBinomial
  real<lower=0> phi_inv; 
}

transformed parameters {
  // Simulation outputs
  real y_sim[n_days, 16];
  
  // Derived metrics
  real yearly_total_RIh_incidence[n_years];
  real annual_new_total_RIh[n_years];
  real annual_new_total_SIh[n_years];
  real total_hai_sih_ratio[n_years];
  real average_total_hai_sih_ratio;

  real community_resistant_proportion[n_years];
  real average_community_resistant_proportion;

  // Pack parameters for ODE solver
  real theta_values[7] = {beta_eh, beta_sc_ec, beta_sc_ac, beta_hh, beta_hc, beta_sh_sih, delta_rcc_ric};

  // 1. Solve ODE
  y_sim = integrate_ode_bdf(sir, y0, t0, ts, theta_values, x_r, x_i);

  // 2. Compute Annual Incidence and Proportions
  for (year in 1:n_years) {
    int start_idx = (year - 1) * 365 + 1;
    int end_idx   = year * 365;

    // Calculate Person-Days at Risk (Susceptible + Colonized)
    vector[365] sh_vec  = to_vector(y_sim[start_idx:end_idx, 1]);
    vector[365] rch_vec = to_vector(y_sim[start_idx:end_idx, 2]);
    real risk_days = sum(sh_vec) + sum(rch_vec) + 1e-6;

    // A. Incidence of Resistant Infections (Hospital)
    // Uses Cumulative Counter (Index 13)
    real prev_cum_RIh    = (year > 1) ? y_sim[start_idx - 1, 13] : y0[13];
    real current_cum_RIh = y_sim[end_idx, 13];
    annual_new_total_RIh[year] = fmax(current_cum_RIh - prev_cum_RIh, 0.0);
    
    // Normalize per 10,000 patient-days
    yearly_total_RIh_incidence[year] = (annual_new_total_RIh[year] / risk_days) * 10000;

    // B. Incidence of Sensitive Infections (Hospital)
    // Uses Cumulative Counter (Index 14)
    real prev_cum_SIh    = (year > 1) ? y_sim[start_idx - 1, 14] : y0[14];
    real current_cum_SIh = y_sim[end_idx, 14];
    annual_new_total_SIh[year] = fmax(current_cum_SIh - prev_cum_SIh, 0.0);

    // C. HAI Ratios
    real total_new_infections = annual_new_total_SIh[year] + annual_new_total_RIh[year];
    if (total_new_infections > 1e-6) {
      total_hai_sih_ratio[year] = annual_new_total_SIh[year] / total_new_infections;
    } else {
      total_hai_sih_ratio[year] = 0;
    }

    // D. Community Metrics
    real prev_cum_RIc    = (year > 1) ? y_sim[start_idx - 1, 15] : y0[15];
    real current_cum_RIc = y_sim[end_idx, 15];
    real new_RIc = fmax(current_cum_RIc - prev_cum_RIc, 0.0);
    
    real prev_cum_SIc    = (year > 1) ? y_sim[start_idx - 1, 16] : y0[16];
    real current_cum_SIc = y_sim[end_idx, 16];
    real new_SIc = fmax(current_cum_SIc - prev_cum_SIc, 0.0);
    
    real total_comm_infections = new_SIc + new_RIc;
    if (total_comm_infections > 1e-6) {
      community_resistant_proportion[year] = new_RIc / total_comm_infections;
    } else {
      community_resistant_proportion[year] = 0;
    }
  }
  
  // Calculate Averages over the simulation period
  average_total_hai_sih_ratio = mean(total_hai_sih_ratio);
  average_community_resistant_proportion = mean(community_resistant_proportion);
}

model {
  // --- Priors ---
  beta_eh       ~ beta(0.19, 0.18);
  beta_sc_ac    ~ beta(0.69, 13.12);
  beta_sc_ec    ~ beta(0.43, 0.35);
  beta_hh       ~ uniform(0.0021, 0.0139);
  beta_hc       ~ uniform(0.0016, 0.011);
  beta_sh_sih   ~ uniform(0, 1);
  delta_rcc_ric ~ uniform(0, 1);
  phi_inv       ~ gamma(50, 500); // Inverse dispersion

  // --- Likelihood (Calibration Targets) ---
  
  // 1. Observed Incidence (Negative Binomial)
  for (i in 1:n_years) {
    real phi = 1.0 / phi_inv;
    RIh_obs[i] ~ neg_binomial_2(fmax(0.0, yearly_total_RIh_incidence[i]), phi);
  }

  // 2. Secondary Targets (Normal approximations)
  average_total_hai_sih_ratio            ~ normal(0.713433333, 0.083628934);
  average_community_resistant_proportion ~ normal(0.090833333, 0.016375795);
}

generated quantities {
  // Export predicted quantities for PPC (Posterior Predictive Checks)
  array[n_years] real pred_yearly_total_RIh_incidence = yearly_total_RIh_incidence;
  array[n_years] real pred_annual_new_total_RIh       = annual_new_total_RIh;
  array[n_years] real pred_annual_new_total_SIh       = annual_new_total_SIh;
  array[n_years] real pred_total_hai_sih_ratio        = total_hai_sih_ratio;
  array[n_years] real pred_community_resistant_proportion = community_resistant_proportion;
  
  real pred_average_total_hai_sih_ratio = average_total_hai_sih_ratio;
  real pred_average_community_resistant_proportion = average_community_resistant_proportion;
}

