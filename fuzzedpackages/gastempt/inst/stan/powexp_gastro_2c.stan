// Fit power exponential gastric emptying curves with Stan
// This is the cholesky variant, which supposedly is more stable
// than variant 2b with corr_matrix instead.
// Covariance between beta and tempt is estimated using LKJ prior.
// However, using random init often fails; use init = 0 instead.
// Initial volume is locally normalized to 1, and denormalized
// in generated quantities for v0.
// Student-t with degrees of freedom that can be set.
// Lower LKJ: lower variance of beta
// Lower values of df: lower variance of beta. Other mainly unaffected

data{
  real lkj; // lkj parameter, see below
  int student_df; // 3 to 9
  int<lower=0> n; // Number of data
  int<lower=0> n_record; // Number of records, used for v0
  int record[n];
  vector[n] minute;
  vector[n] volume;
}

transformed data{
    vector[2] zeros;
    vector[n] volume_1;
    real norm_vol;
    int n_norm;
    zeros[1] = 0;
    zeros[2] = 0;
    n_norm = 0;
    norm_vol = 0;
    // Use mean of initial volume to normalize
    for (i in 1:n){
      if (minute[i] < 5) {
        norm_vol += volume[i];
        n_norm  += 1;
      }
    }
    norm_vol /= n_norm;
    volume_1 = volume/norm_vol;
}


parameters{
  vector <lower=0, upper = 2>[n_record] v0_1;
  vector<lower=0>[2] sigma_record;
  real <lower=0> mu_beta;
  real <lower=0> mu_tempt;
  real<lower=0> sigma;
  cholesky_factor_corr[2] L_rho;
  matrix[2, n_record] z;       // z-scores for constructing varying effects
}

transformed parameters{
  matrix[n_record, 2] cf;
  cf = (diag_pre_multiply(sigma_record, L_rho) * z)';
}

model{
  int  rec;
  real v0r;
  real betar;
  real temptr;
  vector[n] mu;
  // http://www.psychstatistics.com/2014/12/27/d-lkj-priors/
  // Large values, e.g. 70, give priors that prefer low correlations
  // near 0. At 1 flat -1 to 1; lower 1 produces a trough at 0
  L_rho ~ lkj_corr_cholesky(lkj);
  to_vector(z) ~ normal(0,1); // sample z-scores for varying effects

  v0_1  ~ normal(1, 0.3);
  mu_beta ~ normal(1.2, 0.5);
  mu_tempt ~ normal(100, 50);
  sigma_record[1] ~ cauchy(0,20);
  sigma_record[2] ~ cauchy(0,0.4);
  sigma ~ cauchy(0., 0.5);

for (i in 1:n){
   rec = record[i];
   v0r = v0_1[rec];
   temptr = mu_tempt + cf[rec, 1];
   betar = mu_beta + cf[rec, 2];
   mu[i] = v0r*exp(-(minute[i]/temptr)^betar);
  }
  // Using fixed student_t degrees of freedom. Values from 3 (many outliers)
  // to 9 are useful
  volume_1 ~ student_t(student_df, mu, sigma);
}

generated quantities{
  vector[n_record] v0;
  vector[n_record] tempt;
  vector[n_record] beta;
  // Denormalized v0
  v0 = v0_1 * norm_vol;
  for (i in 1:n_record){
    tempt[i] = mu_tempt + cf[i,1];
    beta[i] = mu_beta + cf[i,2];
  }
}
