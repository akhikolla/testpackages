
functions {
real normal_lb_ub_rng(real mu, real sigma, real lb, real ub) {
    real p1 = normal_cdf(lb, mu, sigma);  // cdf with lower bound
    real p2 = normal_cdf(ub, mu, sigma);  // cdf with upper bound
    real u = uniform_rng(p1, p2);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf
}

real get_log_sd(real mu, real sigma){
  real log_sigma = sqrt(log(((sigma^2)/(mu^2)) + 1));
  if (log_sigma <= 0){
    log_sigma = 0.0000000001;
  }
  return(log_sigma);
}

real get_log_mu(real mu, real sigma){
  real log_sigma = sqrt(log(((sigma^2)/(mu^2)) + 1));
  real log_mu = log(mu) - (0.5 * (log_sigma^2));
  return(log_mu);
}


}



data {

                     //means of all parameters

                     real lt_m;
                     real ac_m;
                     real an_m;
                     real ap_m;
                     real Dc_m;
                     real Dn_m;
                     real Dp_m;
                     real linf_m;
                     real k_m;
                     real t0_m;
                     real theta_m;
                     real r_m;
                     real h_m;
                     real lwa_m;
                     real lwb_m;
                     real mdw_m;
                     real v_m;
                     real F0nz_m;
                     real F0pz_m;
                     real Qc_m;
                     real Qn_m;
                     real Qp_m;
                     real alpha_m;
                     real f0_m;

                     //sd of all parameters

                     real lt_sd;
                     real ac_sd;
                     real an_sd;
                     real ap_sd;
                     real Dc_sd;
                     real Dn_sd;
                     real Dp_sd;
                     real linf_sd;
                     real k_sd;
                     real t0_sd;
                     real theta_sd;
                     real r_sd;
                     real h_sd;
                     real lwa_sd;
                     real lwb_sd;
                     real mdw_sd;
                     real v_sd;
                     real F0nz_sd;
                     real F0pz_sd;
                     real Qc_sd;
                     real Qn_sd;
                     real Qp_sd;
                     real alpha_sd;
                     real f0_sd;

                     // correlations ro
                     real ro_Qc_Qn;
                     real ro_Qc_Qp;
                     real ro_Qn_Qp;
                     real ro_Dc_Dn;
                     real ro_Dc_Dp;
                     real ro_Dn_Dp;
                     real ro_lwa_lwb;
                     real ro_alpha_f0;

                     }

transformed data{
  real logQc_m = get_log_mu(Qc_m, Qc_sd);
  real logQc_sd = get_log_sd(Qc_m, Qc_sd);
  real logQn_m = get_log_mu(Qn_m, Qn_sd);
  real logQn_sd = get_log_sd(Qn_m, Qn_sd);
  real logQp_m = get_log_mu(Qp_m, Qp_sd);
  real logQp_sd = get_log_sd(Qp_m, Qp_sd);

  real logDc_m = get_log_mu(Dc_m, Dc_sd);
  real logDc_sd = get_log_sd(Dc_m, Dc_sd);
  real logDn_m = get_log_mu(Dn_m, Dn_sd);
  real logDn_sd = get_log_sd(Dn_m, Dn_sd);
  real logDp_m = get_log_mu(Dp_m, Dp_sd);
  real logDp_sd = get_log_sd(Dp_m, Dp_sd);

  real loglwa_m = get_log_mu(lwa_m, lwa_sd);
  real loglwa_sd = get_log_sd(lwa_m, lwa_sd);
  real loglwb_m = get_log_mu(lwb_m, lwb_sd);
  real loglwb_sd = get_log_sd(lwb_m, lwb_sd);

  real loga_m = get_log_mu(alpha_m, alpha_sd);
  real loga_sd = get_log_sd(alpha_m, alpha_sd);
  real logf0_m = get_log_mu(f0_m, f0_sd);
  real logf0_sd = get_log_sd(f0_m, f0_sd);
}

model{

                     }


generated quantities {

                     //all paramaters

                     real<lower=0.001> lt;
                     real<lower=0> ac;
                     real<lower=0> an;
                     real<lower=0> ap;
                     real<lower=0> Dc;
                     real<lower=0> Dn;
                     real<lower=0> Dp;
                     real<lower=0> linf;
                     real<lower=0> k;
                     real t0;
                     real<lower=0,upper=10> theta;
                     real<lower=0> r;
                     real<lower=0> h;
                     real<lower=0> lwa;
                     real<lower=0> lwb;
                     real<lower=0> mdw;
                     real<lower=0> v;
                     real<lower=0> F0nz;
                     real<lower=0> F0pz;
                     real<lower=0> Qc;
                     real<lower=0> Qn;
                     real<lower=0> Qp;
                     real<lower=0> alpha;
                     real<lower=0> f0;

                     //derived variables

                     real<lower=0> m_max;   //max weight
                     real l1;      //same as lt
                     real a1;      //age at length l1
                     real a2;      // age at next time interval
                     real l2;      // predicted length at age a2
                     real w1;      // weight for l1
                     real w2;      // weight for l2
                     real wd1;     // dry weight equivalent w1
                     real wd2;     // dry weight equivalent w2
                     real Wd;      // gain in dry weight over time step
                     real Ww;      // gain in wet weight over time step

                     real Qc1;      // mass C of fish at l1 in g
                     real Qn1;      // mass N of fish at l1 in g
                     real Qp1;      // mass P of fish at l1 in g
                     real Gc;     // C gain for growth in g
                     real Gn;     // N gain for growth in g
                     real Gp;     // P gain for growth in g

                     real Em;      //cost of growth in J/g
                     real gC_to_J; // conversion factor
                     real Ec;      // combustion energy of biomass (Joules / g)
	                   real Bm;
	                   real B_main;  // maintenance metabolic rate
	                   real B_syn;   // cost of growth
	                   real B_rest;  // resting metabolic rate (Joules / day)
	                   real B_tot;   // assimilation rate in joule per day
	                   real F0c;      // amount of mass C needed for metabolism

                     real F0n;     // N needed for cell renewal
                     real F0p;     // P needed for cell renewal

                     real Sn;
                     real Sp;
                     real Sc;

                    // needed nutrients
                     real st_np;
                     real st_cn;
                     real st_cp;

                    // food
                     real stf_np;
                     real stf_cn;
                     real stf_cp;

                     int lim;      // limiting element


                     real Ic;    // ingestion
                     real In;
                     real Ip;

                     real Wc;
                     real Wn;    // egestion
                     real Wp;

                     real Fn;    // excretion
                     real Fp;

                     real Fc;     // total respiration

                     real Frn;     // leftover excretion
                     real Frp;

                     real IN;      // ingestion in g dry weight
                     real IN_cnp;

          ////////// Body CNP estimation /////////

                     // covariance matrices
                     matrix[3,3] Sigma_Qcnp;
                     matrix[3,3] Sigma_Dcnp;
                     matrix[2,2] Sigma_lw;
                     matrix[2,2] Sigma_ab;
                     // vectors of parameter means
                     vector[3] mu_Qcnp;
                     vector[3] mu_Dcnp;
                     vector[2] mu_lw;
                     vector[2] mu_ab;

                     // vectors of parameter estimates from multinormal sampling
                     vector[3] Qcnp;
                     vector[3] Dcnp;
                     vector[2] lw;
                     vector[2] ab;

                     // fill vectors of parameter means
                     mu_Qcnp[1] = logQc_m;
                     mu_Qcnp[2] = logQn_m;
                     mu_Qcnp[3] = logQp_m;

                     mu_Dcnp[1] = logDc_m;
                     mu_Dcnp[2] = logDn_m;
                     mu_Dcnp[3] = logDp_m;

                     mu_lw[1] = loglwa_m;
                     mu_lw[2] = loglwb_m;

                     mu_ab[1] = loga_m;
                     mu_ab[2] = logf0_m;

                     // construct cov matrices
                     Sigma_Qcnp[1,1] = logQc_sd^2;
                     Sigma_Qcnp[2,2] = logQn_sd^2;
                     Sigma_Qcnp[3,3] = logQp_sd^2;
                     Sigma_Qcnp[1,2] = logQc_sd * logQn_sd * ro_Qc_Qn;
                     Sigma_Qcnp[2,1] = logQc_sd * logQn_sd * ro_Qc_Qn;
                     Sigma_Qcnp[1,3] = logQc_sd * logQp_sd * ro_Qc_Qp;
                     Sigma_Qcnp[3,1] = logQc_sd * logQp_sd * ro_Qc_Qp;
                     Sigma_Qcnp[2,3] = logQn_sd * logQp_sd * ro_Qn_Qp;
                     Sigma_Qcnp[3,2] = logQn_sd * logQp_sd * ro_Qn_Qp;

                     Sigma_Dcnp[1,1] = logDc_sd^2;
                     Sigma_Dcnp[2,2] = logDn_sd^2;
                     Sigma_Dcnp[3,3] = logDp_sd^2;
                     Sigma_Dcnp[1,2] = logDc_sd * logDn_sd * ro_Dc_Dn;
                     Sigma_Dcnp[2,1] = logDc_sd * logDn_sd * ro_Dc_Dn;
                     Sigma_Dcnp[1,3] = logDc_sd * logDp_sd * ro_Dc_Dp;
                     Sigma_Dcnp[3,1] = logDc_sd * logDp_sd * ro_Dc_Dp;
                     Sigma_Dcnp[2,3] = logDn_sd * logDp_sd * ro_Dn_Dp;
                     Sigma_Dcnp[3,2] = logDn_sd * logDp_sd * ro_Dn_Dp;

                     Sigma_lw[1,1] = loglwa_sd^2;
                     Sigma_lw[2,2] = loglwb_sd^2;
                     Sigma_lw[1,2] = loglwa_sd * loglwb_sd * ro_lwa_lwb;
                     Sigma_lw[2,1] = loglwa_sd * loglwb_sd * ro_lwa_lwb;

                     Sigma_ab[1,1] = loga_sd^2;
                     Sigma_ab[2,2] = logf0_sd^2;
                     Sigma_ab[1,2] = loga_sd * logf0_sd * ro_alpha_f0;
                     Sigma_ab[2,1] = loga_sd * logf0_sd * ro_alpha_f0;

                     // sample from multinormal distributions
                     Qcnp = multi_normal_rng(mu_Qcnp, Sigma_Qcnp);
                     Dcnp = multi_normal_rng(mu_Dcnp, Sigma_Dcnp);
                     lw = multi_normal_rng(mu_lw, Sigma_lw);
                     ab = multi_normal_rng(mu_ab, Sigma_ab);

                     // back transform estimates

                     Qc = exp(Qcnp[1]);
                     Qn = exp(Qcnp[2]);
                     Qp = exp(Qcnp[3]);

                     Dc = exp(Dcnp[1]);
                     Dn = exp(Dcnp[2]);
                     Dp = exp(Dcnp[3]);

                     lwa = exp(lw[1]);
                     lwb = exp(lw[2]);

                     alpha = exp(ab[1]);
                     f0 = exp(ab[2]);

                     // Sample other parameters

                      lt = normal_lb_ub_rng(lt_m, lt_sd, 0.5, 1000);
                      ac = normal_lb_ub_rng(ac_m, ac_sd, 0.0001, 1);
                      an = normal_lb_ub_rng(an_m, an_sd, 0.0001, 1);
                      ap = normal_lb_ub_rng(ap_m, ap_sd, 0.0001, 1);
                      linf = normal_lb_ub_rng(linf_m, linf_sd, 1, 1000);
                      k = normal_lb_ub_rng(k_m, k_sd,0.0001,3);
                      t0 = normal_rng(t0_m, t0_sd);
                      theta = normal_lb_ub_rng(theta_m, theta_sd,0.1, 6);
                      r = normal_lb_ub_rng(r_m, r_sd,0.001, 8);
                      h = normal_lb_ub_rng(h_m, h_sd,1, 5);
                      mdw = normal_lb_ub_rng(mdw_m, mdw_sd, 0.001, 1);
                      v = normal_rng(v_m, v_sd);
                      F0nz = normal_lb_ub_rng(F0nz_m, F0nz_sd, 0.000000000000001, 0.1);
                      F0pz = normal_lb_ub_rng(F0pz_m, F0pz_sd, 0.000000000000001, 0.1);

                     //Quantify derived values

                     m_max = lwa * linf^lwb;  //maximum weight based on linf

                     l1  = lt;                                    // lt1, Total length of the fish at the moment
                     w1  = lwa * (l1^lwb);                        // conversion to weight in g
                     wd1 = w1 * mdw;                           // conversion to dry weight in g

                     // Growth per day

                     if (lt < linf){
                     a1  = log(1.0 - (l1/linf))/(-k) + t0;          // Age1, Predicted age of the fish at length lt
                     a2  = a1 + (1.0 / 365);                        // Age2, Age1 + 1 day
                     l2  = linf * (1.0 - exp(-k * (a2 - t0)));      // lt2, Predicted total length at age 2
                     w2  = lwa * (l2^lwb);
                     wd2 = w2 * mdw;
                     Wd  = wd2 - wd1;                             // Growth in dry weight
                     Ww  = w2 - w1;                               // Growth in wet weight
                     }
                     if (lt>= linf){
                     // if bigger than linf, growth is zero
                     a1  = 100;  // arbitrary age, infinity
                     a2  = a1;
                     l2  = l1;
                     w2  = w1;
                     wd2 = wd1;
                     Wd = 0;
                     Ww = 0;
                     }

                     Qc1  = Qc * wd1 / 100;
                     Qn1  = Qn * wd1 / 100;
                     Qp1  = Qp * wd1 / 100;
                     Gc = Qc * Wd / 100;
                     Gn = Qn * Wd / 100;
                     Gp = Qp * Wd / 100;

                     // metabolism

                     	 Em       = exp(4.38 + 0.1032 * log(v) + 0.73 * log(h) + 0.41 * log(r + 1.0));  //cost of growth in J/g

                     	 gC_to_J  = 39e3;                         // conversion factor
                     	 Ec       = 24e3;                         // combustion energy of biomass (Joules / g)
	                     Bm       = f0 * gC_to_J * m_max^(alpha - 1.0);
	                     B_main   = Bm * w1;                       // maintenance metabolic rate
	                     B_syn    = Em * Ww;                       // cost of growth
	                     B_rest   = B_main + B_syn;                // resting metabolic rate (Joules / day)
	                     B_tot    = B_rest * theta;                    // assimilation rate in joule per day
	                     F0c       = B_tot / gC_to_J;               // amount of mass C needed for metabolism

                     // biomass turnover

                       F0n = F0nz * Qn1;                            // N needed for cell renewal
                       F0p = F0pz * Qp1;                            // P needed for cell renewal

                     // Needed ingestion of each element
                       Sn = (Gn + F0n) / an;
                       Sp = (Gp + F0p) / ap;
                       Sc = (Gc + F0c) / ac;

                     // Stoichiometry and determining limiting element
                      // needed nutrients
                       st_np = Sn / Sp;
                       st_cn = Sc / Sn;
                       st_cp = Sc / Sp;

                      // food
                       stf_np = Dn / Dp;
                       stf_cn = Dc / Dn;
                       stf_cp = Dc / Dp;

                     // limiting nutrients: C=1, N=2, P=3

                       if (st_cn > stf_cn && st_cp > stf_cp) {
                                lim = 1;
                       } else if  (st_cn < stf_cn && st_np > stf_np){
                                lim = 2;
                       } else {
                                lim = 3;
                       }

                     // ingestion, based upon limiting element

                       if (lim==3){                       // P limiting
                               Ip = Sp;
                               In = Ip * stf_np;
                               Ic = Ip * stf_cp;
                       } else if (lim==2){                // N limiting
                               In = Sn;
                               Ip = In / stf_np;
                               Ic = In * stf_cn;
                       } else{                            // C limiting
                               Ic = Sc;
                               Ip = Ic / stf_cp;
                               In = Ic / stf_cn;
                       }

                     // egestion
                       Wc = Ic * (1-ac);
                       Wn = In * (1-an);
                       Wp = Ip * (1-ap);

                     // excretion
                       Fn = In - Wn - Gn;
                       Fp = Ip - Wp - Gp;

                     // respiration
                       Fc = Ic - Wc - Gc;

                     // leftover excretion
                       Frn = Fn - F0n;
                       Frp = Fp - F0p;

                     // ingestion in dry weight
                       IN = Ic * 100 / Dc;
                       IN_cnp = Ic + In + Ip;

                     }
