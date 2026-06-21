functions {
  real[] sir(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real Sh = y[1];
    real RCh = y[2];
    real RIh = y[3];
    real SIh = y[4];
    real Sc = y[5];
    real RCc = y[6];
    real RIc = y[7];
    real SIc = y[8];
    real Drh = y[9];
    real Dsh = y[10];
    real Drc = y[11];
    real Dsc = y[12];
    real C_RIh_total = y[13];
    real C_SIh_total = y[14];
    real C_RIc_total = y[15];
    real C_SIc_total = y[16];


    real beta_eh = theta[1];
    real beta_sc_ec = theta[2];
    real beta_sc_ac = theta[3];
    real beta_hh = theta[4];
    real beta_hc = theta[5];
    real beta_sh_sih = theta[6];
    real delta_rcc_ric = theta[7];

    real mu_sih = x_r[1];
    real mu_rih = x_r[2];
    real gamma_sih_sh = x_r[3];
    real gamma_rch_sh = x_r[4];
    real gamma_rih_rch = x_r[5];
    real beta_sc_sic = x_r[6];
    real mu_sic = x_r[7];
    real mu_ric = x_r[8];
    real gamma_sic_sc = x_r[9];
    real gamma_rcc_sc = x_r[10];
    real gamma_ric_rcc = x_r[11];
    real delta_rch_sih = x_r[12];
    real delta_rch_rih = x_r[13];
    real alpha_adm = x_r[14];
    real alpha_dis = x_r[15];
    real mu_b = x_r[16];
    real mu_m = x_r[17];
    real mu_d = x_r[18];
    real delta_rcc_sic = x_r[19];


    real Nh_current = fmax(1e-6, Sh + RCh + RIh + SIh);
    real Nc_current = fmax(1e-6, Sc + RCc + RIc + SIc);
    real N_total_living = Nh_current + Nc_current;

    real dSh_dt =
      - beta_eh * Sh
      - beta_hh * Sh * RCh / Nh_current
      - beta_sh_sih * Sh * SIh / Nh_current
      - alpha_dis * Sh
      - mu_d * Sh
      + gamma_sih_sh * SIh
      + gamma_rch_sh * RCh
      + alpha_adm * Sc;

    real dRCh_dt =
      - delta_rch_rih * RCh
      - gamma_rch_sh * RCh
      - alpha_dis * RCh
      - delta_rch_sih * RCh
      - mu_d * RCh
      + beta_eh * Sh
      + beta_hh * Sh * RCh / Nh_current
      + gamma_rih_rch * RIh
      + alpha_adm * RCc;

    real dRIh_dt =
      - gamma_rih_rch * RIh
      - mu_rih * RIh
      + delta_rch_rih * RCh
      + alpha_adm * RIc;

    real dSIh_dt =
      - gamma_sih_sh * SIh
      - mu_sih * SIh
      + beta_sh_sih * Sh * SIh / Nh_current
      + delta_rch_sih * RCh
      + alpha_adm * SIc;

    real dSc_dt = - beta_sc_ec * Sc - beta_sc_ac * Sc * 0.33 - beta_hc * Sc * RCc / Nc_current - beta_sc_sic * Sc * SIc / Nc_current - alpha_adm * Sc - mu_d * Sc + gamma_sic_sc * SIc + gamma_rcc_sc * RCc + alpha_dis * Sh + mu_b * N_total_living + mu_m * N_total_living;
    real dRCc_dt = -delta_rcc_sic * RCc- delta_rcc_ric * RCc - gamma_rcc_sc * RCc - alpha_adm * RCc - mu_d * RCc + beta_sc_ec * Sc + beta_sc_ac * Sc * 0.33 + beta_hc * Sc * RCc / Nc_current + gamma_ric_rcc * RIc + alpha_dis * RCh;
    real dRIc_dt = - gamma_ric_rcc * RIc - mu_ric * RIc - alpha_adm * RIc + delta_rcc_ric * RCc;
    real dSIc_dt = - gamma_sic_sc * SIc - mu_sic * SIc - alpha_adm * SIc + beta_sc_sic * Sc * SIc / Nc_current+ delta_rcc_sic * RCc;

    real dDrh_dt = mu_rih * RIh;
    real dDsh_dt = mu_sih * SIh;
    real dDrc_dt = mu_ric * RIc;
    real dDsc_dt = mu_sic * SIc;

    real dC_RIh_total_dt = delta_rch_rih * RCh + alpha_adm * RIc;
    real dC_SIh_total_dt = beta_sh_sih * Sh * SIh / Nh_current + delta_rch_sih * RCh + alpha_adm * SIc;
    real dC_RIc_total_dt = delta_rcc_ric * RCc;
    real dC_SIc_total_dt = beta_sc_sic * Sc * SIc / Nc_current+ delta_rcc_sic * RCc;

    

    return {
      dSh_dt, dRCh_dt, dRIh_dt, dSIh_dt,
      dSc_dt, dRCc_dt, dRIc_dt, dSIc_dt,
      dDrh_dt, dDsh_dt, dDrc_dt, dDsc_dt,
      dC_RIh_total_dt, dC_SIh_total_dt, dC_RIc_total_dt, dC_SIc_total_dt
    };
  }
}

data {
  int<lower=1> n_years;
  int<lower=1> n_days;
  real ts[n_days];
  real y0[16];
  real t0;
  int<lower=0> RIh_obs[n_years];

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
  real x_r[19];
  x_r[1] = mu_sih_fixed;
  x_r[2] = mu_rih_fixed;
  x_r[3] = gamma_sih_sh_fixed;
  x_r[4] = gamma_rch_sh_fixed;
  x_r[5] = gamma_rih_rch_fixed;
  x_r[6] = beta_sc_sic_fixed;
  x_r[7] = mu_sic_fixed;
  x_r[8] = mu_ric_fixed;
  x_r[9] = gamma_sic_sc_fixed;
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

  int x_i[0];
}

parameters {
  real<lower=0, upper=1> beta_eh;
  real<lower=0, upper=1> beta_sc_ec;
  real<lower=0, upper=1> beta_sc_ac;
  real<lower=0.0021, upper=0.0139> beta_hh;
  real<lower=0.0016, upper=0.011> beta_hc;
  real<lower=0, upper=1> beta_sh_sih;
  real<lower=0, upper=1> delta_rcc_ric;
  real<lower=0> phi_inv;
}


transformed parameters {
  real y_sim[n_days, 16];

  real yearly_total_RIh_incidence[n_years];
  real annual_new_total_RIh[n_years];
  real annual_new_total_SIh[n_years];
  real total_hai_sih_ratio[n_years];
  real average_total_hai_sih_ratio;

  real community_resistant_proportion[n_years];  // Changed from n_days to n_years
  real average_community_resistant_proportion;

  real theta_values[7] = {beta_eh, beta_sc_ec, beta_sc_ac, beta_hh, beta_hc, beta_sh_sih, delta_rcc_ric};

  y_sim = integrate_ode_bdf(sir, y0, t0, ts, theta_values, x_r, x_i);

for (year in 1:n_years) {
    int start_idx = (year - 1) * 365 + 1;
    int end_idx = year * 365;

    vector[365] sh_vec = to_vector(y_sim[start_idx:end_idx, 1]);
    vector[365] rch_vec = to_vector(y_sim[start_idx:end_idx, 2]);
    real risk_days = sum(sh_vec) + sum(rch_vec) + 1e-6;

    // --- Original Incidence Calculation (for RIh_obs likelihood, uses C_RIh_total: y_sim[,13]) ---
    real prev_cumulative_C_RIh_total = (year > 1) ? y_sim[start_idx - 1, 13] : y0[13];
    real current_cumulative_C_RIh_total = y_sim[end_idx, 13];
    annual_new_total_RIh[year] = fmax(current_cumulative_C_RIh_total - prev_cumulative_C_RIh_total, 0.0);
    yearly_total_RIh_incidence[year] = (annual_new_total_RIh[year] / risk_days) * 10000;

    real prev_cumulative_C_SIh_total = (year > 1) ? y_sim[start_idx - 1, 14] : y0[14];
    real current_cumulative_C_SIh_total = y_sim[end_idx, 14];
    annual_new_total_SIh[year] = fmax(current_cumulative_C_SIh_total - prev_cumulative_C_SIh_total, 0.0);

    // Hospital infection ratio
    real total_new_infections_year = annual_new_total_SIh[year] + annual_new_total_RIh[year];
    if (total_new_infections_year > 1e-6) {
      total_hai_sih_ratio[year] = annual_new_total_SIh[year] / total_new_infections_year;
    } else {
      total_hai_sih_ratio[year] = 0;
    }

    // Community infections
    real prev_cumulative_C_RIc_total = (year > 1) ? y_sim[start_idx - 1, 15] : y0[15];
    real current_cumulative_C_RIc_total = y_sim[end_idx, 15];
    real annual_new_total_RIc = fmax(current_cumulative_C_RIc_total - prev_cumulative_C_RIc_total, 0.0);
    
    real prev_cumulative_C_SIc_total = (year > 1) ? y_sim[start_idx - 1, 16] : y0[16];
    real current_cumulative_C_SIc_total = y_sim[end_idx, 16];
    real annual_new_total_SIc = fmax(current_cumulative_C_SIc_total - prev_cumulative_C_SIc_total, 0.0);
    
    // Community resistant proportion
    real total_new_community_infections_year = annual_new_total_SIc + annual_new_total_RIc;
    if (total_new_community_infections_year > 1e-6) {
      community_resistant_proportion[year] = annual_new_total_RIc / total_new_community_infections_year;
    } else {
      community_resistant_proportion[year] = 0;
    }
  }
  
  average_total_hai_sih_ratio = mean(total_hai_sih_ratio);
  average_community_resistant_proportion = mean(community_resistant_proportion);
}


model {
  beta_eh ~ beta(0.9192455, 0.2944133 );
  beta_sc_ac ~ beta(0.9615812, 24.5673739);
  beta_sc_ec ~ beta(1.0997771, 0.8543970);
  beta_hh ~ uniform(0.0021, 0.0139);
  beta_hc ~ uniform(0.0016, 0.011);
  beta_sh_sih ~ uniform(0, 1);
  delta_rcc_ric ~ uniform (0,1);
  phi_inv ~ gamma(50, 500);


  for (i in 1:n_years) {
    real phi = 1.0 / phi_inv;
    RIh_obs[i] ~ neg_binomial_2(fmax(0.0, yearly_total_RIh_incidence[i]), phi);
  }

  average_total_hai_sih_ratio ~ normal(0.713433333, 0.083628934);
  average_community_resistant_proportion ~ normal(0.090833333, 0.016375795);
}

generated quantities {
  // Predictions - just the model's deterministic outputs
  array[n_years] real pred_yearly_total_RIh_incidence = yearly_total_RIh_incidence;
  array[n_years] real pred_annual_new_total_RIh = annual_new_total_RIh;
  array[n_years] real pred_annual_new_total_SIh = annual_new_total_SIh;
  array[n_years] real pred_total_hai_sih_ratio = total_hai_sih_ratio;
  array[n_years] real pred_community_resistant_proportion = community_resistant_proportion;
  real pred_average_total_hai_sih_ratio = average_total_hai_sih_ratio;
  real pred_average_community_resistant_proportion = average_community_resistant_proportion;
  
  // Optional: Daily time series if needed
  // array[n_days] real pred_Sh = to_array_1d(y_sim[,1]);
  // array[n_days] real pred_RCh = to_array_1d(y_sim[,2]);
  // ... other compartments as needed
}

