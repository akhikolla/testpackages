#include "EMVS.h"


// Functions for both conjugate and independent implementations

double beta(double x, double y) {
  double temp = tgamma(x);
  temp *= tgamma(y);
  temp /= tgamma(x + y);
  return temp;
}

double delogit(double x) {
    x = exp(x);
    x /= (1 + x);
    return x;
}

vec density_norm(vec &x, double mu, double sigma) {
  vec dens = 1/(sigma*sqrt(2*PI)) * ones<vec>(x.n_elem);    
  dens %= exp(-square(x - mu)/(2 * pow(sigma, 2))); 
  return dens;
}


vec density_norm_log(vec &x, double mu, double sigma) {
  vec dens = -0.5*log(sigma*sigma*2*PI) * ones<vec>(x.n_elem);    
  dens += -square(x - mu)/(2 * pow(sigma, 2)); 
  return dens;
}

vec pnorm(vec x) {
  vec value=x/(sqrt(2));
  int n = x.n_elem;
  for (int i = 0; i < n; i++){
    value[i] = erf(value[i]);
    value[i] *= 0.5;
    value[i] += 0.5;
  }
  return value;
}

double erfc_log(double x) {
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0) sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y =  log (((((a5*t + a4)*t) + a3)*t + a2)*t + a1) +log(t) + exp(-x*x);

    return sign*y;
}


vec pnorm_log(vec x) {
  vec value = -x/(sqrt(2));
  int n = x.n_elem;
  for (int i = 0; i < n; i++){
    value[i]=log(0.5) + erfc_log(value[i]);   
  }
  return value;
}

// M-step for p
double M_p(vec &post, double a, double b){
  double res = accu(post) + a - 1;
  res /= (b + a + post.n_elem - 2);
  return res;
}

// Functions for conjugate implementation

mat conj_E_beta_binom(vec &beta_k, double sigma_k, double v0, double v1, double theta, double t) {
  mat stars(beta_k.n_elem, 2);

  // Compute p_star
  vec dens1 = density_norm_log(beta_k, 0.0, sigma_k * sqrt(v1));
  vec dens0 = density_norm_log(beta_k, 0.0, sigma_k * sqrt(v0));

  stars.col(1) = 1/ (1 + pow(1 - theta, t)/pow(theta, t)*exp(t*(dens0 - dens1)));

  // Compute d_star
  stars.col(0) =  stars.col(1)/v1 + (1-stars.col(1))/v0;  
  return stars;
}


vec conj_M_beta(mat &XtY, mat &X, vec &Y, mat &XtX, vec &inv_var) {
  
  if(X.n_cols > X.n_rows){
    mat X_star=X;
    X_star.each_row()/=sqrt(inv_var.t());
    mat Psi=X_star*X_star.t();
    Psi.diag()+=1;
    Psi=solve(Psi,Y);
    X_star.each_row()/=sqrt(inv_var.t());
    return X_star.t()*Psi;
} else {
    mat Psi = XtX;
    Psi.diag()+=inv_var;
    Psi = inv(Psi);
    Psi *= XtY;
    return Psi;
  }
}


double conj_M_sigma(mat &Y, mat &X, vec &beta_k, vec &inv_var, double eta, double lambda){
  vec e = Y - X*beta_k;
  vec e2 = inv_var % beta_k;
  
  double res = as_scalar(e.t() * e);
  res += as_scalar(e2.t() * beta_k);
  res += eta*lambda;
  res /= (X.n_rows + X.n_cols+eta);
  res = sqrt(res);
  return res;
}


double log_prior(uvec& gamma, const string &type, double a, double b, int n){
  double res;
  if(type.compare("fixed") == 0){
    // NOT IMPLEMENTED YET
  } else if (type.compare("betabinomial") == 0){
    double x = gamma.n_elem + a;
    double y = n - gamma.n_elem + b;
    res = log(beta(x,y)) - log(beta(a,b));
    
    // Stirling approximation:
    if(!is_finite(res)){
      res = 0.5*log(2*PI)+(x-0.5)*log(x)+(y-0.5)*log(y)-(x+y-0.5)*log(x+y);
    }
    
  } else if (type.compare("MRF") == 0){
    // NOT IMPLEMENTED YET
  }
  return res;
}


double log_g(uvec &gamma, mat &X, mat &Y, double eta, double lambda, double v0, double v1, const string &type, double a, double b){
  int q = gamma.n_elem;
  double Ssq, log_val;
  int n = Y.n_elem;
  int N = X.n_cols;
  if(q > 0){
    mat X_gamma = X.cols(gamma);
    vec inv_var = ones<vec>(q)/v1;
    mat X_tilde = join_cols(X_gamma, diagmat(sqrt(inv_var)));
    vec XtY= X_gamma.t() * Y;
    mat XtX=X_gamma.t() * X_gamma;
    vec Yresp=Y.col(0);
    vec aux= conj_M_beta(XtY, X_gamma, Yresp, XtX, inv_var);
    double sign=-1;

    Ssq = as_scalar(Y.t()*Y - XtY.t() * aux);
    log_det(log_val,sign,X_tilde.t()*X_tilde);
    log_val*=-0.5;
    log_val += -0.5*q*log(v1);
    log_val -= 0.5*(n+eta)*log(eta*lambda+Ssq);
    log_val += log_prior(gamma,type,a,b,N);
  } else {
    Ssq = as_scalar(Y.t()*Y);
    log_val =  -0.5*(n+eta)*log(eta*lambda+Ssq);
    log_val += log_prior(gamma,type,a,b,N);
  }
  return log_val;
}



// Functions for independent implementation

mat ind_E_beta_binom(vec &beta_k, double v0, double v1, double theta, double t){
  mat stars(beta_k.n_elem, 2);

  // Compute p_star
  vec dens1 = density_norm_log(beta_k, 0.0, sqrt(v1));
  vec dens0 = density_norm_log(beta_k, 0.0, sqrt(v0));

  stars.col(1) = 1/ (1+ pow(1-theta,t)/pow(theta,t)*exp(t*(dens0-dens1)));

  // Compute d_star
  stars.col(0) =  stars.col(1)/v1 + (1-stars.col(1))/v0;  
  return stars;
}


vec ind_M_beta(mat &XtY, mat &X, vec &Y, mat &XtX, vec &inv_var, double sigma_k){
  
  // Multiply sigma by D
  int p = inv_var.n_elem;
  vec sigma_d = ones<vec>(p);
  for (int i = 0; i < p; i++) {
    sigma_d[i] = sigma_k * inv_var[i]; 
  }
  
  if(X.n_cols > X.n_rows){
    mat X_star=X;
    X_star.each_row()/=sqrt(sigma_d.t());
    mat Psi=X_star*X_star.t();
    Psi.diag()+=1;
    Psi=solve(Psi,Y);
    X_star.each_row()/=sqrt(sigma_d.t());
    return X_star.t()*Psi;
} else {
    mat Psi = XtX;
    Psi.diag()+=sigma_d;
    Psi = inv(Psi);
    Psi *= XtY;
    return Psi;
  }
}



double ind_M_sigma(mat &Y, mat &X, vec &beta_k, double eta, double lambda){
  vec e = Y - X*beta_k;

  double res = as_scalar(e.t() * e);
  res += eta*lambda;
  res /= (X.n_rows + eta + 1);
  res = sqrt(res);
  return res;
}

