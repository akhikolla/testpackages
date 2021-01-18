data {
 // Define variables in data
 // Number of observations (an integer)
 int<lower=0> N;
 // Number of fixed parameters
 int<lower=0> p;
  // Number of lambda param
 int<lower=0> q;
 // Number of eps param
 int<lower=0> r;
 // Number of knots i.e. random effects
 int<lower=0> nk;
  // number of columns for each Z matrix
 int zvars[nk>0 ? q+1:0];
 // Variables
 int y[N];
 // fixed effects design matrix
 matrix[N, p] X;
 // random effects design matrix
 matrix[nk > 0 ? N:0, nk>0 ? nk:0] Z;
 // max col of Z
 int max_col;
 // indicator for random intercept;
 int<lower=0,upper=1> randint;
 // indicator for random effects;
 int<lower=0,upper=1> randeff;
 // indicator whether to use QR decomposition
 int<lower=0,upper=1> qr; // 0 = no, 1 = yes
 // indicator whether to split QR decomposition across multiple matrices
 int<lower=0,upper=1> qrsplit;
 // family number:  1=gaussian, 2=binomial, 3=poisson
 int<lower=1, upper=3> famnum;
  // link number
 int linknum;
  // offset
 vector[N] offset;

 ////////////////////////////////////////////////////////
 // hyperpriors
 ////////////////////////////////////////////////////////

 // lambda
 int lambdanum[nk>0 ? (q+1):0];
 int lambda_max_params;
 matrix[nk>0 ? q:0, nk>0 ? lambda_max_params:0] lambda_param;

 // beta
 int betanum[p+1];
 int beta_max_params;
 matrix[p, beta_max_params] beta_param;
}

transformed data {
  // qr for X
  matrix[N, p] Q_x;
  matrix[p, p] R_x;
  matrix[p, p] R_x_inverse;

  // qr for Z
  matrix[N, nk] Q_z;
  matrix[nk, nk] R_z;
  matrix[nk, nk] R_z_inverse;

  // qr for Z with array of matrices
  matrix[N, max_col] qz[q];
  matrix[max_col, max_col] rz[q];
  matrix[max_col, max_col] rzinv[q];


  // thin and scale the QR decomposition
  Q_x = qr_Q(X)[, 1:p] * sqrt(N - 1);
  R_x = qr_R(X)[1:p, ] / sqrt(N - 1);
  R_x_inverse = inverse(R_x);

  // thin and scale the QR decomposition
  if (randeff == 1 || randint == 1) {
    Q_z = qr_Q(Z)[, 1:nk] * sqrt(N - 1);
    R_z = qr_R(Z)[1:nk, ] / sqrt(N - 1);
    R_z_inverse = inverse(R_z);
  }


}

parameters {
 // Define parameters to estimate
 vector[p] theta_b;
 vector[nk] tau;
 vector<lower=0>[q] lambda;
}

transformed parameters {
  vector[nk] theta_u;
  vector[p] beta;
  vector[nk] u;

  if (randeff == 1 || randint == 1) {
    for (l4 in 1:1) {
      int i = 1;
      for (j4 in 1:q) {
        for (k4 in 1:zvars[j4]) {
          theta_u[i] = tau[i] * lambda[j4];
          i = i + 1;
        }
      }
    }
  }
  
  if (qr == 1) {
    beta = R_x_inverse * theta_b; // coefficients on x
    if (randeff == 1 || randint == 1) {
      u = R_z_inverse * theta_u;      
    }
  } else {
    beta = theta_b;
    if (randeff == 1 || randint == 1) {
      u = theta_u;      
    }
  }

}

model {
 // Prior
 if (randint == 1 || randeff == 1) {
    for (k3 in 1:q) {
     if (lambdanum[k3] == 1) {
       lambda[k3] ~ normal(lambda_param[k3, 1], lambda_param[k3, 2]);
     } else if (lambdanum[k3] == 2) {
       lambda[k3] ~ student_t(lambda_param[k3, 1], lambda_param[k3, 2], lambda_param[k3, 3]);
     }
   }
 }


   for (k1 in 1:p) {
     if (betanum[k1] == 1) {
       theta_b[k1] ~ normal(beta_param[k1, 1], beta_param[k1, 2]);
     } else if (betanum[k1] == 2) {
       theta_b[k1] ~ student_t(beta_param[k1, 1], beta_param[k1, 2], beta_param[k1, 3]);
     }
   }

   if (randeff == 1 || randint == 1) {
      for (k2 in 1:nk)
       tau[k2] ~ normal(0, 1);    
   }

 // Binomial
 if (randeff == 1 || randint == 1) {
   if (famnum == 2) {
     // logit link
     if (linknum == 4) {
        if (qr == 1) {
          y ~ bernoulli_logit(Q_x*theta_b + Q_z*theta_u + offset);
        } else {
          y ~ bernoulli_logit(X*theta_b + Z*theta_u + offset);
        }
      // probit link
     } else if (linknum == 5) {
       if (qr == 1) {
          y ~ bernoulli(Phi(Q_x*theta_b + Q_z*theta_u + offset));
        } else {
          y ~ bernoulli(Phi(X*theta_b + Z*theta_u + offset));
        }
      // cauchit link
     } else if (linknum == 6) {
       if (qr == 1) {
          y ~ bernoulli(atan(Q_x*theta_b + Q_z*theta_u + offset) / pi() + 0.5);
        } else {
          y ~ bernoulli(atan(X*theta_b + Z*theta_u + offset) / pi() + 0.5);
        }
      // log link
     } else if (linknum == 2) {
       if (qr == 1) {
          y ~ bernoulli(exp(Q_x*theta_b + Q_z*theta_u + offset));
        } else {
          y ~ bernoulli(exp(X*theta_b + Z*theta_u + offset));
        }
      // cloglog link
     } else if (linknum == 7) {
       if (qr == 1) {
          y ~ bernoulli(inv_cloglog(Q_x*theta_b + Q_z*theta_u + offset));
        } else {
          y ~ bernoulli(inv_cloglog(X*theta_b + Z*theta_u + offset));
        }
     }
    // Poisson
   } else if (famnum == 3) {
     // log link
     if (linknum == 2) {
        if (qr == 1) {
          y ~ poisson_log(Q_x*theta_b + Q_z*theta_u + offset);
        } else {
          y ~ poisson_log(X*theta_b + Z*theta_u + offset);
        }
      // identity link
     } else if (linknum == 1) {
       if (qr == 1) {
          y ~ poisson(Q_x*theta_b + Q_z*theta_u + offset);
        } else {
          y ~ poisson(X*theta_b + Z*theta_u + offset);
        }
      // sqrt link
     } else if (linknum == 8) {
       if (qr == 1) {
          y ~ poisson(square(Q_x*theta_b + Q_z*theta_u + offset));
        } else {
          y ~ poisson(square(X*theta_b + Z*theta_u + offset));
        }
     }
  
   }
 } else {
   if (famnum == 2) {
     // logit link
     if (linknum == 4) {
        if (qr == 1) {
          y ~ bernoulli_logit(Q_x*theta_b + offset);
        } else {
          y ~ bernoulli_logit(X*theta_b + offset);
        }
      // probit link
     } else if (linknum == 5) {
       if (qr == 1) {
          y ~ bernoulli(Phi(Q_x*theta_b + offset));
        } else {
          y ~ bernoulli(Phi(X*theta_b + offset));
        }
      // cauchit link
     } else if (linknum == 6) {
       if (qr == 1) {
          y ~ bernoulli(atan(Q_x*theta_b + offset) / pi() + 0.5);
        } else {
          y ~ bernoulli(atan(X*theta_b + offset) / pi() + 0.5);
        }
      // log link
     } else if (linknum == 2) {
       if (qr == 1) {
          y ~ bernoulli(exp(Q_x*theta_b + offset));
        } else {
          y ~ bernoulli(exp(X*theta_b + offset));
        }
      // cloglog link
     } else if (linknum == 7) {
       if (qr == 1) {
          y ~ bernoulli(inv_cloglog(Q_x*theta_b + offset));
        } else {
          y ~ bernoulli(inv_cloglog(X*theta_b + offset));
        }
     }
    // Poisson
   } else if (famnum == 3) {
     // log link
     if (linknum == 2) {
        if (qr == 1) {
          y ~ poisson_log(Q_x*theta_b + offset);
        } else {
          y ~ poisson_log(X*theta_b + offset);
        }
      // identity link
     } else if (linknum == 1) {
       if (qr == 1) {
          y ~ poisson(Q_x*theta_b + offset);
        } else {
          y ~ poisson(X*theta_b + offset);
        }
      // sqrt link
     } else if (linknum == 8) {
       if (qr == 1) {
          y ~ poisson(square(Q_x*theta_b + offset));
        } else {
          y ~ poisson(square(X*theta_b + offset));
        }
     }
  
   }
 }

}

generated quantities {
  vector[N] log_lik;
  
  // extract log_lik
  if (randint == 1 || randeff == 1) {
   for (n in 1:N) {
     // Binomial
     if (famnum == 2) {
       // logit link
       if (linknum == 4) {
          log_lik[n] = bernoulli_logit_lpmf(y[n] | X[n, ]*beta + Z[n,]*u);
        // probit link
       } else if (linknum == 5) {
          log_lik[n] = bernoulli_lpmf( y[n] | Phi(X[n, ]*beta + Z[n,]*u));
        // cauchit link
       } else if (linknum == 6) {
         log_lik[n] = bernoulli_lpmf( y[n] | atan(X[n, ]*beta + Z[n,]*u)/pi() + 0.5);
        // log link
       } else if (linknum == 2) {
         log_lik[n] = bernoulli_lpmf( y[n] | exp(X[n, ]*beta + Z[n,]*u));
        // cloglog link
       } else if (linknum == 7) {
         log_lik[n] = bernoulli_lpmf( y[n] | inv_cloglog(X[n, ]*beta + Z[n,]*u));
       }
    
    
      // Poisson
     } else if (famnum == 3) {
       // log link
       if (linknum == 2) {
         log_lik[n] = poisson_log_lpmf(y[n] | X[n, ]*beta + Z[n,]*u);
        // identity link
       } else if (linknum == 1) {
         log_lik[n] = poisson_lpmf(y[n] | X[n, ]*beta + Z[n,]*u);
        // sqrt link
       } else if (linknum == 8) {
         log_lik[n] = poisson_lpmf(y[n] | square(X[n, ]*beta + Z[n,]*u));
       }
     }
   }
  } else {
  for (n in 1:N) {
     // Binomial
     if (famnum == 2) {
       // logit link
       if (linknum == 4) {
          log_lik[n] = bernoulli_logit_lpmf(y[n] | X[n, ]*beta);
        // probit link
       } else if (linknum == 5) {
          log_lik[n] = bernoulli_lpmf( y[n] | Phi(X[n, ]*beta));
        // cauchit link
       } else if (linknum == 6) {
         log_lik[n] = bernoulli_lpmf( y[n] | atan(X[n, ]*beta)/pi() + 0.5);
        // log link
       } else if (linknum == 2) {
         log_lik[n] = bernoulli_lpmf( y[n] | exp(X[n, ]*beta));
        // cloglog link
       } else if (linknum == 7) {
         log_lik[n] = bernoulli_lpmf( y[n] | inv_cloglog(X[n, ]*beta));
       }
    
    
      // Poisson
     } else if (famnum == 3) {
       // log link
       if (linknum == 2) {
         log_lik[n] = poisson_log_lpmf(y[n] | X[n, ]*beta);
        // identity link
       } else if (linknum == 1) {
         log_lik[n] = poisson_lpmf(y[n] | X[n, ]*beta);
        // sqrt link
       } else if (linknum == 8) {
         log_lik[n] = poisson_lpmf(y[n] | square(X[n, ]*beta));
       }
     }
   }
  }
}

