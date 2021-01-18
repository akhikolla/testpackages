functions {
  real normal_lb_ub_rng(real mu, real sigma, real lb, real ub) {
    real p1 = normal_cdf(lb, mu, sigma);  // cdf with lower bound
    real p2 = normal_cdf(ub, mu, sigma);  // cdf with upper bound
    real u = uniform_rng(p1, p2);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf
}

real get_B0_rng(real troph_m, real troph_sd, real temp, real b0_m, real b0_sd) {
    real troph = normal_lb_ub_rng(troph_m, troph_sd,0,5);
    real b0 = normal_lb_ub_rng(b0_m, b0_sd,0,1);
    real Er     = 0.70; // activation energy (eV)
		real Ei     = 2.43; // inactivation parameter (eV)
		real tl     = 0.81; // trophic level slope
		real lnb0tc = log(b0); // metabolic normalization independent of mass, temperature, and trophic level (g C / troph^0.81 / g^0.77 / day)
		real k      = 8.62e-5; // Boltzmann constant (eV / K)
		real Tc     = 293.15; // arbitrary absolute temperature (K)
		real Tk     = temp + 273.15; // temperature (K)
		real Topt   = 306.4; // optimum temperature (K)
		real kt     = (1 / (k * Tc)) - (1 / (k * Tk)); // (1 / eV)
		real ktop   = (1 / (k * Topt)) - (1 / (k * Tk)); // (1 / eV)
    real lnB0  = lnb0tc + tl * log(troph) + (Er * kt) - log(1 + ((Er / (Ei - Er)) * exp(Ei * ktop)));
    real B0 = exp(lnB0);
    return(B0);
    }
}

data {
real troph_m;
real troph_sd;
real temp;
real b0_m;
real b0_sd;
}

model { }

generated quantities {
real B0;
B0 = get_B0_rng(troph_m,  troph_sd,  temp,  b0_m,  b0_sd);
}



