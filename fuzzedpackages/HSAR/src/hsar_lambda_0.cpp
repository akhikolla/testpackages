#include <iostream>
#include <math.h>   
#include <RcppArmadillo.h>
#include "diagnostics.h"
#include "gibbs_method.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// Estimate the y hat 
mat y_hat_hsar_lambda_0(const mat& X, const mat& betas, double rho, const sp_mat& W, const sp_mat& Z,const mat& Mus  ){
  
  mat Xb = X*betas.t();
  int n = X.n_rows;
  
  sp_mat I_sp = speye<sp_mat>(n,n);
  
  //mat invA = spsolve(A,I_nn, "lapack");
  sp_mat invA = I_sp+rho*W+pow(rho,2)*(W*W)+pow(rho,3)*(W*W*W);
  
  mat Zg = Z* Mus.t();
  
  mat yhat = invA*(Xb+Zg);

  return yhat;
}

// Log likelihood function of HSAR
double HSAR_loglikelihood_lambda_0(const mat& X, const mat& y, double rho, 
                                const mat& betas, const mat& us, 
                                const vec& Unum, int Utotal, double sigma2e, 
                                const mat& detval, const sp_mat& W){
  
  int n = X.n_rows;
  
  //find log( det(A) )
  uvec indt = find(detval.col(0) > rho, 1, "first");
  double logdetA = detval(indt[0],1);
  
  sp_mat I_sp = speye<sp_mat>(n,n);
  sp_mat A = I_sp - rho * W;
    
  mat Ay = A*y;
  
  mat Xb = X*trans(betas);
  
  mat Zu;
  mat nus = trans(us);
  for(int j=0;j<Utotal;j++){
    mat temp_add = repmat(nus.row(j),Unum[j], 1 );
      
    Zu.insert_rows( Zu.n_rows, temp_add);
  }
 
  mat crosprod_AymXbmZu= trans(Ay-Xb-Zu)*(Ay-Xb-Zu);
    
  double log_lik_sar = (-n/2)*(log(2*datum::pi)+log( pow(sigma2e,2))) 
  + logdetA - crosprod_AymXbmZu(0,0)/(2*pow(sigma2e,2));
  
  return log_lik_sar;
}

// [[Rcpp::export]]

List hsar_cpp_arma_lambda_0( arma::mat X, arma::vec y, arma::sp_mat W, 
                      arma::sp_mat Z, arma::mat detval, arma::vec Unum, 
                      int burnin, int Nsim, int thinning,
                      float rho_start, float sigma2e_start, float sigma2u_start
                               , arma::vec betas_start){
  
  //Starting the MCMC SET UP
  //arma_rng::set_seed(124);
  
  int n = X.n_rows;
  int p = X.n_cols;
  
  int Utotal (Unum.n_elem );
  //Prior distribution specifications
  //For Betas
  vec M0 = betas_start;
  
  mat T0 = 100.0 * eye<mat>(p,p);
  
  //completely non-informative priors
  int c0(0.01);
  int d0(0.01);
  int a0(0.01);
  int b0(0.01);
  
  //Store MCMC results
  mat Betas = zeros(Nsim-burnin, p);
  mat Us = zeros(Nsim-burnin, Utotal);
  
  vec sigma2e = zeros(Nsim-burnin);
  vec sigma2u = zeros(Nsim-burnin);
  
  vec rho = zeros(Nsim-burnin);
  
  //initial values for model parameters
  float sigma2e_i(sigma2e_start), sigma2e_ip1(sigma2e_start);
  float sigma2u_i(sigma2u_start), sigma2u_ip1(sigma2u_start);
  
  float rho_i(rho_start), rho_ip1(rho_start);
  
  int ce( n/2 + c0 );
  int au( Utotal/2 + a0 );
  
  //Fixed matrix manipulations during the MCMC loops
  
  mat XTX = trans(X) * X;
  mat invT0 = inv(T0);
  
  mat T0M0 = invT0 * M0;
  mat tX = trans(X);
  mat Zfull(Z);
  mat tZ = trans(Zfull);
  
  //some fixed values when updating rho
  mat beta0 = solve(X,y);
  
  mat e0 = y-X*beta0;
  mat e0e0 = trans(e0) * e0;
  
  vec Wy = W *y;
  mat betad = solve(X,Wy);
  mat ed = Wy - X * betad;
  
  mat eded = trans(ed) *ed;
  mat e0ed = trans(e0) *ed;
  
  // initialise A
  sp_mat I_sp_A = speye<sp_mat>(n,n);
  sp_mat A = I_sp_A - rho_i * W;

  // MCMC updating  
  mat us = zeros(Utotal,1);
  for(int i=1;i<=Nsim*thinning;i++) {
    // Gibbs sampler for updating Betas
    mat VV = (1.0/sigma2e_i) *XTX + invT0;
    // We use Cholesky decomption to inverse the covariance matrix
    mat vBetas = inv_sympd(VV) ;//chol2inv(chol(VV))
    
    // Define A=I-rho*W 
    A = I_sp_A - rho_i * W;
    
    vec Ay = A*y;
    vec ZUs = Z*us ;
    
    mat mBetas = vBetas * (tX * (1.0/sigma2e_i)*(Ay-ZUs)+T0M0);    
    
    // When the number of independent variables is large, drawing from multivariate 
    // distribution could be time-comsuming
    mat cholV = trans(chol(vBetas));
    // draw form p-dimensitional independent norm distribution
    mat betas = Rcpp::rnorm(p);
    betas = mBetas + cholV * betas;
  
    // update U. us is spatially independent.

    vec vU = (sigma2e_i * sigma2u_i) /(sigma2u_i*Unum + sigma2e_i); 
     
    vec Xb = X*betas;
    
    vec mU = (1.0/sigma2e_i)*vU%(trans(Z)*(Ay - Xb));
  
    us = zeros(Utotal,1);
     
    for(int j=0;j<Utotal;j++){
      us(j,0) = Rcpp::rnorm(1,mU[j],sqrt(vU[j]))[0];
    }
    
    // Gibbs sampler for updating sigma2e
    mat Zu;
    for(int j=0;j<Utotal;j++){
      mat temp_add = repmat(us.row(j),Unum[j], 1 );
      
      Zu.insert_rows( Zu.n_rows, temp_add);
    }
    
    mat e = Ay - Zu -Xb;
    mat de = 0.5 * trans(e) * e + d0;
    
    sigma2e_ip1 = 1/Rf_rgamma(ce,1/de(0,0));
    
   // Gibbs sampler for updating sigma2u
   mat bu = 0.5 * trans(us) * us + b0;
   sigma2u_ip1 = 1/Rf_rgamma(au,1/bu(0,0));
   
   // Giddy Gibbs integration and inverse sampling for rho
    mat betau = solve(X,Zu);
    mat eu = Zu-X*betau;
    
    mat eueu = trans(eu) * eu;
    mat e0eu = trans(e0) * eu;
    mat edeu = trans(ed) * eu;
    
    rho_ip1 = HSAR_draw_rho(detval, e0e0, eded, eueu, e0ed, e0eu, edeu, sigma2e_ip1);  
    
    rho_i = rho_ip1;
    sigma2e_i = sigma2e_ip1;sigma2u_i = sigma2u_ip1;
    
    if(i>(burnin*thinning)) 
    {
      if (i%thinning==0) 
      {
        Betas.row(int(i-burnin*thinning)/thinning-1) = trans(betas);
        Us.row(int(i-burnin*thinning)/thinning-1) = trans(us);
        
        sigma2e[int(i-burnin*thinning)/thinning-1] = sigma2e_i;
        sigma2u[int(i-burnin*thinning)/thinning-1] = sigma2u_i;
        
        rho[int(i-burnin*thinning)/thinning-1] = rho_i;
      }
    }
    
  }
  
  // Diagnostics
  vec log_lik_samples = zeros(Nsim-burnin);
  
  for(int i=0;i<(Nsim-burnin);i++)
  {
  log_lik_samples[i] = HSAR_loglikelihood_lambda_0( X, y, rho[i], Betas.row(i), 
                                    Us.row(i),Unum,Utotal,
                                    sigma2e[i], detval, W );
  }
  
  double log_lik_mean_theta = HSAR_loglikelihood_lambda_0( X, y, mean( rho ),
                              mean( Betas ), 
                              mean( Us ), Unum,Utotal,
                             mean( sigma2e ), detval, 
                             W);
                             
  double dic, pd;
  diagnostic_dic_pd(log_lik_samples,log_lik_mean_theta, dic, pd);
  
  double r2 = diagnostic_Rsquared(y, y_hat_hsar_lambda_0(X, mean( Betas ), mean( rho ), W, Z ,mean( Us ) ));
  
  mat direct, indirect, total;
  diagnostic_impacts( mean( Betas ), mean( rho ),W , direct, indirect, total);
          
  return List ::create( Named("cbetas")= Betas,
                             Named("Mbetas")= mean( Betas ), 
                             Named("SDbetas") = stddev( Betas ),
                             Named("Mrho")= mean( rho ), 
                             Named("SDrho") = stddev( rho ),
                             Named("Msigma2e")= mean( sigma2e ), 
                             Named("SDsigma2e") = stddev( sigma2e ),
                             Named("Msigma2u")= mean( sigma2u ), 
                             Named("SDsigma2u") = stddev( sigma2u ),
                             Named("Mus")= mean( Us ), 
                             Named("SDus") = stddev( Us ),
                             Named("DIC") = dic,
                             Named("pD") = pd,
                             Named("Log_Likelihood") = log_lik_mean_theta,
                             Named("R_Squared") = r2,
                             Named("impact_direct") = direct,
                             Named("impact_idirect") = indirect,
                             Named("impact_total") = total
                             ); 
    
}


