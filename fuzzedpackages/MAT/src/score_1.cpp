#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

// [[Rcpp::export]]

List SCORE_cpp(arma::mat ipar, 
                      arma::mat resp_full,  
                      const int p,
                      arma::mat sigma, 
                      const int maxIter = 30,
                      const double conv = 0.001,
                      const double D = 1.7,
                      bool Fisher = true) { 
  
  //score input
  
  int ni = ipar.n_rows;
  int nj = resp_full.n_rows;
  arma::mat TH = arma::mat(nj,p);
  arma::mat SE = arma::mat(nj,p);
  
  arma::colvec pos_ths = arma::colvec(p); arma::colvec pre_ths;
  int iter ;
  bool converged ;
  arma::rowvec resp ; 
  
  arma::mat delta ; arma::mat abs_delta ;
  arma::vec d2_diag ; arma::mat d2_inv ;
  bool conv_status ;
  
  //dLL input
  
  arma::mat sigma_inv = arma::mat(p,p);     
  arma::rowvec dll = arma::zeros<arma::rowvec>(p); arma::rowvec deriv1 = arma::rowvec(p);
  arma::mat a = ipar.cols(0,p-1);
  arma::colvec d = ipar.col(p);
  arma::colvec c = ipar.col(p+1);
  arma::mat P;
  double num;
  double u;
  arma::rowvec w;
  bool Bayesian = true ;
  
  //makeFI input
  
  arma::mat FI = arma::zeros<arma::mat>(p,p); arma::mat deriv2;
    arma::mat FI_temp = arma::zeros<arma::mat>(p,p);
    arma::mat cf;
    bool addsigma = true;
    
  //makeHessian input
  
  arma::mat H = arma::zeros<arma::mat>(p,p);
    arma::mat H_temp = arma::zeros<arma::mat>(p,p);
  
  for (int j=0; j<nj; j++) { 
  
  pos_ths = arma::zeros<arma::colvec>(p);
  iter = 0 ;
  converged = false ;
  resp = resp_full.row(j);
  
  while ((iter < maxIter) & (converged == false)) {
    
    iter++ ;
    pre_ths = pos_ths ;
    
    {
    
    // calculate first derivative
      
    dll = arma::zeros<arma::rowvec>(p);
    
    if (p==1) { sigma_inv=1; }
    else { sigma_inv = inv(sigma); }
    
    for (int i=0; i<ni; i++) {   
    u = arma::as_scalar(resp.col(i));
    if ((u==1) | (u==0))
      {
      P = c.row(i) + (1-c.row(i))/(1+exp(-D*(a.row(i)*pre_ths + d.row(i))));
      num = arma::as_scalar((P-c.row(i))*(u-P)/((1-c.row(i))*P));
      dll = dll + a.row(i)*num;
      }
    }
    
    if (Bayesian == true)
    {
      for (int h=0; h<p; h++) {
      w = arma::zeros<arma::rowvec>(p);
      w.col(h) = 1;
      dll.col(h) = dll.col(h) - w*sigma_inv*pre_ths;
      }
    }
    
    dll = arma::as_scalar(D)*dll ;
    
    deriv1 = dll ;
                      }
    
    //calculate second derivative
    
    if (Fisher == true) {
    
    FI = arma::zeros<arma::mat>(p,p);
    FI_temp = arma::zeros<arma::mat>(p,p);
    
    for (int i=0; i<ni; i++) {
    
    P = c.row(i) + (1-c.row(i))/(1+exp(-D*(a.row(i)*pre_ths + d.row(i))));
    cf = (1-P)*pow(P-c.row(i),2.0)/(P*pow(1-c.row(i),2.0));
    FI_temp = trans(a.row(i))*a.row(i);
    num = arma::as_scalar(pow(D,2.0)*cf);
    FI_temp = num*FI_temp;
    FI = FI + FI_temp;
    
    }
    
    if (addsigma == true)
    {
    FI = FI + inv(sigma);
    }
    
    deriv2 = -FI ;
                      }
      
    else {
      
      H = arma::zeros<arma::mat>(p,p);
      H_temp = arma::zeros<arma::mat>(p,p);
      
      for (int i=0; i<ni; i++) {
    u = arma::as_scalar(resp.col(i));
    if ((u==1) | (u==0))
      {    
      P = c.row(i) + (1-c.row(i))/(1+exp(-D*(a.row(i)*pre_ths + d.row(i))));
      cf = arma::as_scalar((1-P)*(P-c.row(i))*(c.row(i)*u-pow(P,2.0))/(pow(P,2.0)*pow(1-c.row(i),2.0)));
      H_temp = trans(a.row(i))*a.row(i);
      H_temp = arma::as_scalar(pow(D,2.0)*cf)*H_temp;
      H = H + H_temp;  
      }
    }
    if (addsigma == true)
    {
    H = H - inv(sigma);
    }
    
    deriv2 = H;
    }
    
    //calculate delta
    
    delta = inv(deriv2)*trans(deriv1);
    pos_ths = pre_ths - delta ;
    abs_delta = abs(delta);
    conv_status = all(vectorise(abs_delta) < conv );
    if (conv_status == true) { converged = true; }
    
                      }
    
    //makeHessian for deriv2
    
    H = arma::zeros<arma::mat>(p,p);
      H_temp = arma::zeros<arma::mat>(p,p);
      
      for (int i=0; i<ni; i++) {
    u = arma::as_scalar(resp.col(i));
    if ((u==1) | (u==0))
      {    
      P = c.row(i) + (1-c.row(i))/(1+exp(-D*(a.row(i)*pos_ths + d.row(i))));
      cf = arma::as_scalar((1-P)*(P-c.row(i))*(c.row(i)*u-pow(P,2.0))/(pow(P,2.0)*pow(1-c.row(i),2.0)));
      H_temp = trans(a.row(i))*a.row(i);
      H_temp = arma::as_scalar(pow(D,2.0)*cf)*H_temp;
      H = H + H_temp;  
      }
    }
    if (addsigma == true)
    {
    H = H - inv(sigma);
    }
    
    deriv2 = H;
    
    TH.row(j) = trans(pos_ths) ;
    d2_inv = inv(deriv2);
    d2_diag = d2_inv.diag();
    SE.row(j) = trans(pow(abs(d2_diag),0.5));
    
  }
    
    return List::create(Named("theta")=TH,Named("SE")=SE);
    
  }
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
