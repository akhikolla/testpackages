// distrib.cpp

#include "dmbc.h"

// Probability function for the product of independent bernoulli random variables
void dprodber(double* prob, const int* d, const double* pi, int m, int logscale){
  double prob_tmp = 0;
  if(logscale){
    *prob = 0;
    for(int i = 0; i < m; i++){
      prob_tmp = d[i]*log(pi[i]) + (1 - d[i])*log(1 - pi[i]);
      if(ISNAN(prob_tmp)){
        prob_tmp = log(pow(pi[i], d[i])*pow(1 - pi[i], 1 - d[i]));
      }
      *prob += prob_tmp;
    }
  }
  else{
    *prob = 1;  
    for(int i = 0; i < m; i++){
      prob_tmp = pow(pi[i], d[i])*pow(1 - pi[i], 1 - d[i]);
      *prob *= prob_tmp;
    }
  }
}

// Density function for a multivariate normal random variable
void dmultinorm(double* dens, const double* x, const double* mean, const double* sigma, int n, int p, int logscale){
  arma::mat x_arma(n, p);
  arma::vec mean_arma(p);
  arma::mat sigma_arma(p, p);
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      x_arma(i, j) = x[n*j + i];
    }
    mean_arma(j) = mean[j];
    for(int k = 0; k < p; k++){
      sigma_arma(j, k) = sigma[p*k + j];
    }
  }

  double logdet, sign;
  arma::log_det(logdet, sign, sigma_arma);
  arma::vec distval = mahalanobis(x_arma, mean_arma, sigma_arma);
  for(int i = 0; i < n; i++){
    dens[i] = -(p * log(2*M_PI) + logdet + distval[i])/2.0;
  }

  if(logscale == 0){
    for(int i = 0; i < n; i++){
      dens[i] = exp(dens[i]);
    }
  }
}

// Multivariate normal random deviates
void rmultinorm(double* dev, int n, const double* mean, const double* sigma, int p){
  arma::vec mean_arma(p);
  arma::mat sigma_arma(p, p);
  for(int j = 0; j < p; j++){
    mean_arma(j) = mean[j];
    for(int i = 0; i < p; i++){
      sigma_arma(i, j) = sigma[i + p*j];
    }
  }

  arma::mat V(p, p);
  arma::mat U(p, p);
  arma::vec d(p);
  arma::svd(U, d, V, sigma_arma);
  arma::mat U_trans = trans(U);

  arma::mat tmp(p, p);
  for(int i = 0; i < p; i++){
    for(int j = 0; j < p; j++){
      tmp(i, j) = U_trans(i, j)*sqrt(d(i));
    }
  }
  arma::mat sigmaSVDsqrt = trans(V * tmp);
  arma::mat normdev(n, p);
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      normdev(i, j) = R::rnorm(0, 1);
    }
  }
  arma::mat out_tmp = normdev * sigmaSVDsqrt;
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      dev[n*j + i] = out_tmp(i, j) + mean_arma(j);
    }
  }
}

// Density function for an inverse gamma random variable
void dinvgamma(double* dens, const double* x, const double alpha, const double beta, int n, int logscale){
  if((alpha <= 0) || (beta <= 0)){
    error("alpha (shape) and beta (scale) parameters in dinvgamma() need to be both strictly positive.\n");
  }

  double lbeta = log(beta);
  double lgalpha = R::lgammafn(alpha);
  for(int i = 0; i < n; i++){
    dens[i] = alpha * lbeta - lgalpha - (alpha + 1) * log(x[i]) - (beta/x[i]);
  }

  if(logscale == 0){
    for(int i = 0; i < n; i++){
      dens[i] = exp(dens[i]);
    }
  }
}

// Inverse gamma random deviates
void rinvgamma(double* dev, int n, const double alpha, const double beta){
  if((alpha <= 0) || (beta <= 0)){
    error("alpha (shape) and beta (scale) parameters in rinvgamma() need to be both strictly positive.\n");
  }

  for(int i = 0; i < n; i++){
    dev[i] = 1.0/R::rgamma(alpha, 1.0/beta);
  }
}

// Density function for a Dirichlet random variable
void ddirichlet(double* dens, const double* x, const double* par, int n, int p, int logscale){
  double tmp = 0;
  arma::mat x_arma(n, p);
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      x_arma(i, j) = x[n*j + i];
    }
  }

  tmp = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      tmp += x_arma(i, j);
      if(x_arma(i, j) < 0 || x_arma(i, j) > 1){
        error("some elements of x outside the [0, 1] range.\n");
      }
    }
    if(std::fabs(tmp - 1.0) >= sqrt(std::numeric_limits<double>::epsilon())){
      error("some rows of x sum to a value different from 1.\n");
    }
  }

  double logD = 0, salpha = 0, s = 0;
  for(int j = 0; j < p; j++){
    logD += R::lgammafn(par[j]);
    salpha += par[j];
  }
  logD = logD - R::lgammafn(salpha);
  for(int i = 0; i < n; i++){
    s = 0;
    for(int j = 0; j < p; j++){
      s += (par[j] - 1)*log(x_arma(i, j));
    }
    dens[i] = s - logD;
    if(logscale == 0){
      dens[i] = exp(dens[i]);
    }
  }
}

// Dirichlet random deviates
void rdirichlet(double* dev, int n, const double* par, int p){
  arma::mat z(n, p, arma::fill::zeros);
  arma::vec s(n, arma::fill::zeros);

  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      z(i, j) = R::rgamma(par[j], 1.0);
      s(i) += z(i, j);
    }
  }
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      dev[n*j + i] = z(i, j)/s(i);
    }
  }
}

// Loglikelihood of the model (contribution to the loglikelihood from a single cluster)
// (d is a vector which stacks the dissimilarities by individual)
void loglik_rbmds_binom(double* loglik, const int* d, const double* z, const double alpha, int n, int p, int S){
  int m = n*(n - 1)/2;

  double* delta = new double[m];
  dist(delta, z, n, p);

  int* diss = new int[m];
  for(int k = 0; k < m; k++){
    diss[k] = 0;
    for(int s = 0; s < S; s++){
      diss[k] += d[k + m*s];
    }
    delta[k] = alpha + delta[k];
  }

  double* pi = new double[m];
  expit(pi, delta, m);

  *loglik = 0;
  for(int k = 0; k < m; k++){
    *loglik += log(pow(pi[k], diss[k])*pow(1 - pi[k], S - diss[k]));
  }

  delete[] delta;
  delete[] diss;
  delete[] pi;
}

// Overall loglikelihood of the model
// (d is a vector which stacks the dissimilarities by individual)
void loglik_dmbc(double* loglik, const int* d, const double* z, const double* alpha, const double* sigma2, const double* lambda, const int* x, int n, int p, int S, int G, const char* family){
  int m = n*(n - 1)/2;
  int ind = 0;

  int* ng = new int[G];
  for(int g = 0; g < G; g++){
    ng[g] = 0;
  }
  tableC(ng, x, S, G);

  double* z_g = new double[n*p];
  int* D_g = new int[m*S];
  double* log_f_g = new double[G];
  *loglik = 0;
  for(int g = 0; g < G; g++){
    for(int j = 0; j < p; j++){
      for(int i = 0; i < n; i++){
        z_g[i + n*j] = z[i + n*j + n*p*g];
      }
    }
    for(int s = 0; s < S; s++){
      if(x[s] == (g + 1)){
        for(int k = 0; k < m; k++){
          D_g[ind] = d[k + m*s];
          ind++;
        }
      }
    }
    if(std::strcmp(family, "binomial") == 0){
      loglik_rbmds_binom(&log_f_g[g], D_g, z_g, alpha[g], n, p, ng[g]);
    }
    *loglik += ng[g]*log(lambda[g]) + log_f_g[g];
    ind = 0;
  }

  delete[] ng;
  delete[] z_g;
  delete[] D_g;
  delete[] log_f_g;
}
