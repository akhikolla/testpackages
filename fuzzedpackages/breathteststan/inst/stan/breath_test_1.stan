// Simple non-hierarchial ungrouped fit of Breath Test Curves to Exponential beta

data{
  int<lower=0> n; // Number of data
  int<lower=0> n_record; // Number of records
  int<lower=1> student_t_df; // using Gaussian for student_t_df >= 10
  real<lower=0> dose;
  int<lower=0> pat_group_i[n];
  vector<lower=0>[n] minute;
  vector<lower=-30>[n] pdr;
}

parameters{
  vector[n_record] m_raw;
  real<lower=0> mu_m;
  real<lower=0> sigma_m;

  vector[n_record] k_raw;
  real<lower=0> mu_k;
  real<lower=0> sigma_k;

  vector[n_record] beta_raw;
  real<lower=0> mu_beta;
  real<lower=0> sigma_beta;

  real <lower=0> sigma;
}


transformed parameters{
  vector<lower=0>[n_record] m;
  vector<lower=0>[n_record] k;
  vector<lower=0>[n_record] beta;

  // Re-parametrization
  m    = mu_m + sigma_m * m_raw;
  k    = mu_k + sigma_k * k_raw;
  beta = mu_beta + sigma_beta * beta_raw;
}

model {
  // Note: the x_raw parameters all have normal(0,1) here
  m_raw ~ normal(0, 1);
  mu_m ~ normal(40,30);
  sigma_m ~ cauchy(0,10);

  k_raw ~ normal(0, 1);
  mu_k ~ lognormal(-5, 2);
  sigma_k ~ lognormal(-7, 2);

  beta_raw ~ normal(0, 1);
  mu_beta ~ normal(2,0.5);
  sigma_beta ~ cauchy(0,2);

  sigma ~ cauchy(0,5);
  { // Dummy block to hide pdr1[n]
    vector[n] pdr1;
    for (i in 1:n){
      int rec;
      real mn;
      real exp_ktn;
      real kn;
      real betan;
      rec = pat_group_i[i];
      mn  =  m[rec];
      kn = k[rec];
      exp_ktn = exp(-kn* minute[i]);
      betan = beta[rec];
      pdr1[i] = dose*mn*kn*betan*exp_ktn * pow(1 - exp_ktn,(betan -1));
    }
    if (student_t_df < 10 )
      pdr ~ student_t(student_t_df, pdr1, sigma);
    else
      pdr ~ normal(pdr1, sigma);
  }
}
