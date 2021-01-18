#include "gibbs_method.h"

// draw rho-values using griddy Gibbs and inversion
double HSAR_draw_rho(const mat& detval,const mat& e0e0,const mat& eded, const mat& eueu,const mat& e0ed,const mat& e0eu, const mat& edeu, double sig){
  
  double rho_draw = 0;
  int nrho = detval.n_rows;
  
  vec rho_exist = detval.col(0);
  vec log_detrho = detval.col(1);
  
  vec iota = ones<vec>(nrho);
  
  // Calculate Log-S(e(rho*)) given both rho* and u*
  
  vec S_rho = e0e0(0,0)*iota + pow(rho_exist, 2) *eded(0,0) +eueu(0,0)*iota - 2*rho_exist*e0ed(0,0) - 2*e0eu(0,0)*iota +2*rho_exist*edeu;
  
  // Calculate the Log-density
  vec log_den = log_detrho - (1.0/(2.0*sig)) * S_rho; 
  
  double adj = log_den.max();
  vec t = rep(adj, nrho);
  log_den = log_den - t ;
  
  // the density
  vec den = exp(log_den);
  // the interval separating rho is h=0.001
  double h = 0.001;
  
  // Integration to calculate the normalized constant
  // using the  trapezoid rule
  double ISUM = h* ( sum(den)-den[0]/2.0-den[nrho-1]/2.0  );
  vec norm_den = (1.0/ISUM)*den;
  
  // cumulative density
  vec cumu_den = cumsum(norm_den);
  
  // Inverse sampling
  double rnd = runif(1,0,1)(0)*sum(norm_den);
  uvec ind = find(cumu_den <= rnd) ;
  int idraw = max(ind);
  
  if(idraw >= 0 && idraw < nrho) 
    rho_draw = rho_exist(idraw);
  
  return rho_draw;
}

// draw lambda-values using griddy Gibbs and inversion, see Smith and LeSage (2004)
double HSAR_draw_lambda(const mat& detvalM,const mat& uu,const mat& uMu, const mat& uMMu, double sig){
  
  double lambda_draw = 0;
  int nlambda = detvalM.n_rows;
  
  vec lambda_exist = detvalM.col(0);
  vec log_detlambda = detvalM.col(1);
  
  vec iota = ones<vec>(nlambda);
  
  // Calculate Log-S(e(rho*)) given both rho* and u*
  vec S_lambda = uu(0,0)*iota + - 2*lambda_exist*uMu + pow(lambda_exist, 2) *uMMu(0,0);
  
  // Calculate the Log-density
  vec log_den = log_detlambda - (1.0/(2.0*sig)) * S_lambda; 
  double adj = log_den.max();
  vec t = rep(adj, nlambda);
  log_den = log_den - t ;
  
  // the density
  vec den = exp(log_den);
  // the interval separating rho is h=0.001
  double h = 0.001;
  
  // Integration to calculate the normalized constant
  // using the  trapezoid rule
  double ISUM = h* ( sum(den)-den[0]/2.0-den[nlambda-1]/2.0  );
  vec norm_den = (1.0/ISUM)*den;
  
  // cumulative density
  vec cumu_den = cumsum(norm_den);
  
  // Inverse sampling
  double rnd = runif(1,0,1)(0)*sum(norm_den);
  uvec ind = find(cumu_den <= rnd) ;
  int idraw = max(ind);
  
  if(idraw >= 0 && idraw < nlambda) 
    lambda_draw = lambda_exist(idraw);
  
  return lambda_draw;
}

