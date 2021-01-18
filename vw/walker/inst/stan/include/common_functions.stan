// Functions for Kalman filter and smoother for dynamic regression model
// note that these functions are not fully optimised yet

// univariate Kalman filter for RW1+RW2 model, returns the log-likelihood
real gaussian_filter(vector y, int[] y_miss, vector a1, matrix P1, real Ht, 
  matrix Tt, matrix Rt, matrix xreg, vector gamma2_y) {
  
  int k = rows(xreg);
  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;

  vector[m] x = a1;
  matrix[m, m] P = P1;

  for (t in 1:n) {
    real F = quad_form(P[1:k, 1:k], xreg[, t]) + gamma2_y[t] * Ht;
    
    if (y_miss[t] == 0 && F > 1.0e-12) { // protect against numerical issues
      real v = y[t] - dot_product(xreg[, t], head(x, k));
      vector[m] K = P[1:m, 1:k] * xreg[, t] / F;
      x = Tt * (x + K * v);
      P = quad_form_sym(P - K * K' * F, Tt');
      for (i in 1:m) {
         P[i, i] += Rt[i, t];
      }
      loglik -= 0.5 * (log(F) + v * v / F);
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt');
      for (i in 1:m) {
         P[i, i] += Rt[i, t];
      }
    }
  }
  return loglik;
  
}

matrix gaussian_smoother(vector y, int[] y_miss, vector a1, matrix P1, real Ht, 
  matrix Tt, matrix Rt, matrix xreg,vector gamma2_y) {

  int k = rows(xreg);
  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;
  vector[m] x = a1;
  matrix[m, m] P = P1;
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n + 1] r;
  vector[m] tmpr;
  
  for (t in 1:n) {
    
    F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + gamma2_y[t] * Ht;
    
    if (y_miss[t] == 0 && F[t] > 1.0e-12) {
      v[t] = y[t] - dot_product(xreg[, t], head(x, k));
      K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
      x = Tt * (x + K[,t] * v[t]);
      P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt');
      for (i in 1:m) {
         P[i, i] += Rt[i, t];
      }
      loglik -= 0.5 * (log(F[t]) + v[t] * v[t] / F[t]);
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt');
      for (i in 1:m) {
         P[i, i] += Rt[i, t];
      }
    }
  }

  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[, t+1];
    if(y_miss[t] == 0 && F[t] > 1.0e-12) {
      vector[m] tmp2 = rep_vector(0.0, m);
      tmp2[1:k] = xreg[, t];
      r[ ,t] = tmp2 * v[t] / F[t] + (Tt - Tt * K[,t] * tmp2')' * tmp;
    } else {
      r[,t] = Tt' * tmp;
    }
  }
  
  tmpr = r[,1];
  r[,1] = a1 + P1 * tmpr;
  for (t in 2:n) {
    vector[m] tmp = r[,t-1];
    vector[m] tmp2 = r[,t];
    r[,t] = Tt * tmp + Rt[, t] .* tmp2;
  }
  return r[1:m, 1:n];
}



// univariate Kalman filter & smoother for non-gaussian model, 
// returns the log-likelihood of the corresponding approximating Gaussian model
// and a extra correction term
vector glm_approx_loglik(vector y, int[] y_miss, vector a1, matrix P1, vector Ht, 
  matrix Tt, matrix Rt, matrix xreg, int distribution, int[] u, 
  vector y_original, vector xbeta_fixed) {

  int k = rows(xreg);
  int n = rows(y);
  int m = rows(a1);
  vector[2] loglik = rep_vector(0.0, 2);
  vector[m] x = a1;
  matrix[m, m] P = P1;
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n+1] r;
  vector[m] tmpr;
  
  for (t in 1:n) {
    
    F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + Ht[t];
    
    if (y_miss[t] == 0 && F[t] > 1.0e-12) {
      v[t] = y[t] - dot_product(xreg[, t], head(x, k));
      K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
      x = Tt * (x + K[,t] * v[t]);
      P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt');
      for (i in 1:m) {
         P[i, i] += Rt[i, t];
      }
      loglik[1] -= 0.5 * (log(F[t]) + v[t] * v[t] / F[t]);
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt');
      for (i in 1:m) {
         P[i, i] += Rt[i, t];
      }
    }
  }

  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[, t+1];
    if(y_miss[t] == 0 && F[t] > 1.0e-12) {
      vector[m] tmp2 = rep_vector(0.0, m);
      tmp2[1:k] = xreg[, t];
      r[ ,t] = tmp2 * v[t] / F[t] + (Tt - Tt * K[,t] * tmp2')' * tmp;
    } else {
      r[,t] = Tt' * tmp;
    }
  }
  
  tmpr = r[,1];
  r[,1] = a1 + P1 * tmpr;
  for (t in 2:n) {
    vector[m] tmp = r[,t-1];
    vector[m] tmp2 = r[,t];
    r[,t] = Tt * tmp + Rt[, t] .* tmp2;
  }

  // add a correction term
  if (distribution == 1) {
    for(t in 1:n) {
      if (y_miss[t] == 0){
        real xbeta_rw = dot_product(xreg[,t], r[1:k, t]);
        loglik[2] += y_original[t] * (xbeta_rw + xbeta_fixed[t]) - 
        u[t] * exp(xbeta_rw + xbeta_fixed[t]) +
          0.5 * (y[t] - xbeta_rw)^2 / Ht[t];
      }
    }
  } else {
    for(t in 1:n) {
      if (y_miss[t] == 0){
        real xbeta_rw = dot_product(xreg[,t], r[1:k, t]);
        loglik[2] += y_original[t] * (xbeta_rw + xbeta_fixed[t]) - 
        u[t] * log1p(exp(xbeta_rw + xbeta_fixed[t])) +
          0.5 * (y[t] - xbeta_rw)^2 / Ht[t];
      }
    }
  }
  return loglik;
}


// univariate Kalman filter & smoother for non-gaussian model, 
// returns the log-likelihood of the corresponding approximating Gaussian model
// and a extra correction term
matrix glm_approx_smoother(vector y, int[] y_miss, vector a1, matrix P1, vector Ht, 
  matrix Tt, matrix Rt, matrix xreg) {

  int k = rows(xreg);
  int n = rows(y);
  int m = rows(a1);
  vector[m] x = a1;
  matrix[m, m] P = P1;
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n+1] r;
  vector[m] tmpr;
  
  for (t in 1:n) {
    
    F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + Ht[t];
    
    if (y_miss[t] == 0 && F[t] > 1.0e-12) {
      v[t] = y[t] - dot_product(xreg[, t], head(x, k));
      K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
      x = Tt * (x + K[,t] * v[t]);
      P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt');
      for (i in 1:m) {
         P[i, i] += Rt[i, t];
      }
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt');
      for (i in 1:m) {
         P[i, i] += Rt[i, t];
      }
    }
  }

  r[, n + 1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[, t+1];
    if(y_miss[t] == 0 && F[t] > 1.0e-12) {
      vector[m] tmp2 = rep_vector(0.0, m);
      tmp2[1:k] = xreg[, t];
      r[, t] = tmp2 * v[t] / F[t] + (Tt - Tt * K[,t] * tmp2')' * tmp;
    } else {
      r[, t] = Tt' * tmp;
    }
  }
  
  tmpr = r[,1];
  r[,1] = a1 + P1 * tmpr;
  for (t in 2:n) {
    vector[m] tmp = r[,t-1];
    vector[m] tmp2 = r[,t];
    r[,t] = Tt * tmp + Rt[, t] .* tmp2;
  }

  return r[1:m, 1:n];
}

