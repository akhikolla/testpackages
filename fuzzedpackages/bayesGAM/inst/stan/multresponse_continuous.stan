// Continuous multresponse
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
 // Number of random effect parameters or knots
 int<lower=0> nk;
 int<lower=0> ny;
 matrix[N, ny] y;
 // Variables
 // fixed effects design matrix
 matrix[N, p] X;
 // max col of Z
 int max_col;
 // number of columns for each Z matrix
 int zvars[nk>0 ? q+1:0];
 // indicator for random intercept;
 int<lower=0,upper=1> randint;
 // indicator for random effects;
 int<lower=0,upper=1> randeff;
 // number of columns for Z intercept
 int<lower=0> nrandint;
 // number of columns for Z nonparametric
 int<lower=0> nnp;
 // random effects design matrix
 matrix[randeff ? N:0, randeff ? nnp:0] Znp;
 // random effects random intercept matrix
 matrix[randint ? N:0, randint ? nrandint:0] Zint;
 // matrix[N, max_col] Zarray[q];
 matrix[nk>0 ? N:0, nk>0 ? nk:0] Z;
 // indicator whether to use QR decomposition
 int<lower=0,upper=1> qr; // 0 = no, 1 = yes
 // indicator whether to split QR decomposition across multiple matrices
 int<lower=0,upper=1> qrsplit;
 // indicator of multivariate independence
 int<lower=0,upper=1> mvindep;
 // family number:  1=gaussian, 2=binomial, 3=poisson
 int<lower=1, upper=3> famnum;
 // link number
 int linknum;
 // offset - TODO currently unused
 // vector[N] offset;

////////////////////////////////////////////////////////
 // hyperpriors
 ////////////////////////////////////////////////////////

 // epsilon
 int epsnum[r+1];
 int eps_max_params;
 matrix[r, eps_max_params] eps_param;

 // beta
 int betanum[p*r+1];
 int beta_max_params;
 matrix[p*r, beta_max_params] beta_param;

  // lambda for nonparametric
 int lambdanum[nk>0 ? q*r+1:0];
 int lambda_max_params;
 matrix[nk>0 ? q*r:0, nk>0 ? lambda_max_params:0] lambda_param;

 // number of off-diagonal
 int a_num_offdiagonal;
 int anum[randint ? (a_num_offdiagonal+1):0];
 int a_max_params;
 matrix[a_num_offdiagonal, randint ? a_max_params:0] a_param;
}


transformed data {
  // qr for X
  matrix[N, p] Q_x;
  matrix[p, p] R_x;
  matrix[p, p] R_x_inverse;

  // qr for Z knots matrices
  matrix[N, nnp] Q_z;
  matrix[nnp, nnp] R_z;
  matrix[nnp, nnp] R_z_inverse;

  int q_reff;
  int q_rint;

  // thin and scale the QR decomposition X
  Q_x = qr_Q(X)[, 1:p] * sqrt(N - 1);
  R_x = qr_R(X)[1:p, ] / sqrt(N - 1);
  R_x_inverse = inverse(R_x);

  // thin and scale the QR decomposition
  if (randeff == 1) {
    Q_z = qr_Q(Znp)[, 1:nnp] * sqrt(N - 1);
    R_z = qr_R(Znp)[1:nnp, ] / sqrt(N - 1);
    R_z_inverse = inverse(R_z);    
  } 
  
  if (randint == 1) {
    q_rint = 1;
    q_reff = q-1;
  } else {
    q_rint = 0;
    q_reff = q;
  }


}

parameters {
 // Define parameters to estimate
 vector[p] theta_b[ny];

 matrix[randint ? nrandint:0, randint ? ny:0] rint_u_transpose;
 vector<lower=0>[randint == 1 ? ny:0] lambda_rint;

 // nonparametric, if any
 vector[randeff ? nnp:0] tau[ny];

 // residual sd
 vector<lower=0>[randeff == 1 ? q_reff:0] lambda_reff[ny];
 vector<lower=0>[r] eps;

 // a parameters
 vector[mvindep ? 0 : a_num_offdiagonal] a;
}

transformed parameters {
  vector[randeff ? nnp:0] theta_u[ny];
  vector[p] beta[ny];
  vector[nk>0 ? nrandint+nnp:0] u[ny];
  vector[randeff ? nnp:0] reff_u[ny];
  vector[randint ? nrandint:0] rint_u[ny];
  
  /////////////////////////////////////////////////////////////
  // random intercept
  // create diagonal covariance matrix for now
  matrix[randint ? ny:0, randint ? ny:0] sigma_u_random;

  // local block
  if (randint == 1) {
     for (ll2 in 1:1)
      {
        matrix[ny, ny] L;
        matrix[ny, ny] Dhalf;
    
        // assign LDLT decomposition
        Dhalf = diag_matrix(lambda_rint);
        L = diag_matrix(rep_vector(1.0, ny));
        for (ll in 1:1) {
          int iter = 1;
           for (ii in 1:ny) {
            for (jj in 1:ny) {
              if (jj > ii) {
                if (mvindep == 1) {
                  L[jj, ii] = 0;
                } else {
                  L[jj, ii] = a[iter];  
                }
                iter = iter + 1;
              }
            }
          }
        }
        
      // sigma_u_random = L * Dhalf * Dhalf * (L');
      sigma_u_random = tcrossprod(L * Dhalf);
      }
 
  }
  

  /////////////////////////////////////////////////////////////
  // alternate approach
  if (randeff == 1) {
    int zindex = 0;
    if (randint == 1) {
      zindex = 1;
    }
    for (l4 in 1:ny) {
      int i = 1;
      for (j4 in 1:q_reff) {
        for (k4 in 1:zvars[j4+zindex]) {
          theta_u[l4][i] = tau[l4][i] * lambda_reff[l4][j4];
          i = i + 1;
        }
      }
    }    
  }

  if (qr == 1) {
      for (jj in 1:ny) {
        beta[jj] = R_x_inverse * theta_b[jj];
        if (randeff == 1) {
           reff_u[jj] = R_z_inverse * theta_u[jj];
        }

      }
    } else {
      for (jj in 1:ny) {
        beta[jj] = theta_b[jj];
       if (randeff == 1) {
         reff_u[jj] = theta_u[jj];
       }
      }
  }


  for (jj in 1:ny) {
    if (randint == 1) {
      for (kk in 1:nrandint) {
       u[jj][kk] = rint_u_transpose[kk][jj];
       rint_u[jj][kk] = rint_u_transpose[kk][jj];
      }      
    }
    if (randeff == 1) {
      for (ll in 1:nnp) {
        u[jj][nrandint+ll] = reff_u[jj][ll]; 
      }      
    }
  }
  

}

model {
 // off diagonal w prior input
 if (mvindep == 0) {
   for (jj in 1:a_num_offdiagonal) {
    if (anum[jj] == 1) {
       a[jj] ~ normal(a_param[jj, 1], a_param[jj, 2]);
     } else if (anum[jj] == 2) {
       a[jj] ~ student_t(a_param[jj, 1], a_param[jj, 2], a_param[jj, 3]);
     }
   }
 }

 for (jj in 1:nrandint) {
    rint_u_transpose[jj] ~ multi_normal(rep_vector(0.0, ny), sigma_u_random);
 }

 for (j1 in 1:r) {
   for (k1 in 1:p) {
     if (betanum[k1*j1] == 1) {
       theta_b[j1, k1] ~ normal(beta_param[k1*j1, 1], beta_param[k1*j1, 2]);
     } else if (betanum[k1] == 2) {
       theta_b[j1, k1] ~ student_t(beta_param[k1*j1, 1], beta_param[k1*j1, 2], beta_param[k1*j1, 3]);
     }
   }
  }

  // nested loop for multvariate response
  for (j1 in 1:r) {
    if (randint == 1) {
      if (lambdanum[1] == 1) {
        lambda_rint[j1] ~ normal(lambda_param[j1, 1], lambda_param[j1, 2]);
      } else if (lambdanum[1] == 2) {
        lambda_rint[j1] ~ student_t(lambda_param[j1, 1], lambda_param[j1, 2], lambda_param[j1, 3]);
      }      
    }

    if (randeff == 1) {
      for (k1 in 1:q_reff) {
       if (lambdanum[k1*j1] == 1) {
         lambda_reff[j1, k1] ~ normal(lambda_param[k1*j1, 1], lambda_param[k1*j1, 2]);
       } else if (lambdanum[k1*j1] == 2) {
         lambda_reff[j1, k1] ~ student_t(lambda_param[k1*j1, 1], lambda_param[k1*j1, 2], lambda_param[k1*j1, 3]);
       }
     }

    }
  }
  
 for (jj in 1:ny)
   tau[jj] ~ normal(0, 1);

  // nested loop for multvariate response
  if (famnum == 1) {

      vector[N] yhat[ny];
      // multivariate response
      for (jj in 1:ny) {
        if (randint == 1 && randeff == 1) {
         yhat[jj] = Q_x*theta_b[jj] +  Q_z[jj]*theta_u[jj] + Zint*col(rint_u_transpose, jj);
        } else if (randint == 0 && randeff == 1) {
         yhat[jj] = Q_x*theta_b[jj] +  Q_z[jj]*theta_u[jj];
        } else if (randint == 1 && randeff == 0) {
          yhat[jj] = Q_x*theta_b[jj] +  Zint*col(rint_u_transpose, jj);
        } else if (randint == 0 && randeff == 0) {
          yhat[jj] = Q_x*theta_b[jj];
        }
      }
  
     for (k2 in 1:r) {
        if (epsnum[k2] == 1) {
         eps[k2] ~ normal(eps_param[k2, 1], eps_param[k2, 2]);
       } else if (epsnum[k2] == 2) {
         eps[k2] ~ student_t(eps_param[k2, 1], eps_param[k2, 2], eps_param[k2, 3]);
       }
     }

       // identity link
      for (jj in 1:ny) {
         if (linknum == 1) {
           y[, jj] ~ normal(yhat[jj], eps[jj]);
         }
         // log link
         else if (linknum == 2) {
           y[, jj] ~ normal(exp(yhat[jj]), eps[jj]);
         }
         // inverse link
         else if (linknum == 3) {
           y[, jj] ~ normal(inv(yhat[jj]), eps[jj]);
         }
      }
    }
}

generated quantities {
  vector[randint ? ny:0] dhalf_inv;
  matrix[randint ? ny:0, randint ? ny:0] sigma_u_correlation;
  vector[N] log_lik[ny];

  // correlation matrix
  if (randint == 1) {
    dhalf_inv = diagonal(sigma_u_random); 
    for (jj in 1:ny) {
      dhalf_inv[jj] = 1 / sqrt(dhalf_inv[jj]);
    }
  
    sigma_u_correlation = quad_form_diag(sigma_u_random, dhalf_inv);
  }

  // extract log lik
  if (randint > 0 || randeff > 0) {
    for (jj in 1:ny) {
      for (n in 1:N) {
           if (linknum == 1) {
             log_lik[jj][n] = normal_lpdf(y[n] | X[n, ]*beta[jj] + Z[n,]*u[jj], eps[jj]);
           }
           // log link
           else if (linknum == 2) {
             log_lik[jj][n] = normal_lpdf(exp(y[n]) | X[n, ]*beta[jj] + Z[n,]*u[jj], eps[jj]);
           }
           // inverse link
           else if (linknum == 3) {
             log_lik[jj][n] = normal_lpdf(inv(y[n]) | X[n, ]*beta[jj] + Z[n,]*u[jj], eps[jj]);
           }
      }
    }    
  } else {
    for (jj in 1:ny) {
      for (n in 1:N) {
           if (linknum == 1) {
             log_lik[jj][n] = normal_lpdf(y[n] | X[n, ]*beta[jj], eps[jj]);
           }
           // log link
           else if (linknum == 2) {
             log_lik[jj][n] = normal_lpdf(exp(y[n]) | X[n, ]*beta[jj], eps[jj]);
           }
           // inverse link
           else if (linknum == 3) {
             log_lik[jj][n] = normal_lpdf(inv(y[n]) | X[n, ]*beta[jj], eps[jj]);
           }
      }
    }    
  }

}
