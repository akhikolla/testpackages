/*******************************R-callable function***************************************/
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 
#include <iostream>
#include <RcppArmadillo.h>
#include "spm_integral_md.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

/*============GLOBAL VARIABLES===========*/
/*=====END OF GLOBAL VARIABLES===========*/

/*============FUNCTION DEFINITIONS*/



long double mu_quad_linear(double t, 
                      arma::mat y1, 
                      arma::mat gamma1, 
                      arma::mat fH, 
                      arma::mat f1H, 
                      double mu0H, 
                      double thetaH, 
                      arma::mat QH, arma::mat Q1, bool gomp) {
  arma::mat hfH;
  arma::mat hf1H;
  long double mu0Ht;
  arma::mat mu;
  
  hfH = fH - y1;
  hf1H = f1H - y1;
  if(gomp) {
    mu0Ht = mu0H*exp(thetaH*t);
  } else {
    mu0Ht = mu0H;
  }
  arma::mat QH_gamma1 = QH*gamma1;
  //mu = mu0Ht + (hfH.t()*QH)*hfH + arma::sum((QH_gamma1).diag());
  //mu = mu0Ht + hfH.t()*Q1 + Q1.t()*hfH + arma::sum((QH_gamma1).diag());
  mu = mu0Ht + (hfH.t()*QH)*hfH + arma::sum((QH_gamma1).diag()) + hfH.t()*Q1 + Q1.t()*hfH;
  
  return (long double)mu(0,0);
}

long double mu_linear(double t, 
                           arma::mat y1, 
                           arma::mat gamma1, 
                           arma::mat fH, 
                           arma::mat f1H, 
                           double mu0H, 
                           double thetaH, 
                           arma::mat QH, bool gomp) {
  arma::mat hfH;
  arma::mat hf1H;
  long double mu0Ht;
  arma::mat mu;
  
  hfH = fH - y1;
  hf1H = f1H - y1;
  if(gomp) {
    mu0Ht = mu0H*exp(thetaH*t);
  } else {
    mu0Ht = mu0H;
  }
  arma::mat QH_gamma1 = QH*gamma1;
  //mu = mu0Ht + hfH.t()*Q1 + Q1.t()*hfH + arma::sum((QH_gamma1).diag());
  
  mu = mu0Ht + arma::sum((QH_gamma1).diag()) + hfH.t()*QH + QH.t()*hfH;
  
  return (long double)mu(0,0);
}


//Calculating m (y[1]) & gamma(y[2]):
void func1_quad_linear(arma::mat *res, double t, arma::mat *y, arma::mat fH, arma::mat f1H, arma::mat aH, 
                  arma::mat bH, arma::mat QH, double theta) {
  
  arma::mat hfH, hf1H, dy1, dy2;
  //hfH = fH.t() - y[0]; 
  //hf1H = f1H.t() - y[0];
  hfH = fH - y[0]; 
  hf1H = f1H - y[0];
  
  res[0] = -1.00 * (aH*hf1H) + 2.00 * ((y[1]*Q(t, QH, theta))*hfH);
  res[1] = aH*y[1] + y[1]*aH.t() + bH*bH.t() - 2.00 * ((y[1]*Q(t, QH, theta))*y[1]);
  //std::cout << "theta: " << theta << " QH:" << QH << " y[0]:" << y[0] << " y[1]:" << y[1] << " Q:" << Q(t, QH, theta) << " res0:" << res[0] << " res1:" << res[1]<< std::endl;
}



RcppExport SEXP complikMD_quadratic_linear(SEXP dat, 
                                 SEXP n, 
                                 SEXP m, 
                                 SEXP ah, 
                                 SEXP f1h, 
                                 SEXP qh, 
                                 SEXP bh, 
                                 SEXP fh, 
                                 SEXP mu0h, 
                                 SEXP thetah, 
                                 SEXP q1h,
                                 SEXP k, 
                                 SEXP pinv_tol, 
                                 SEXP gomp_) {
  
  long N = as<long>(n); // Number of rows
  long M = as<long>(m); // Number of columns
  
  arma::mat aH = as<arma::mat>(ah);
  arma::mat f1H = as<arma::mat>(f1h);
  arma::mat QH = as<arma::mat>(qh);
  arma::mat Q1 = as<arma::mat>(q1h);
  arma::mat bH = as<arma::mat>(bh);
  arma::mat fH = as<arma::mat>(fh);
  double mu0H = as<double>(mu0h); 
  double thetaH  = as<double>(thetah);
  int dim = as<int>(k);
  bool gomp = as<bool>(gomp_);
  
  //Actual data set
  arma::mat dd = as<arma::mat>(dat);  
  double ptol = as<double>(pinv_tol);
  double L; // Likelihood
  
  //End of data loading
  arma::mat *out, *k1ar, *yfin, *ytmp, *k2ar, *k3ar, *k4ar;
  out = new arma::mat[2];
  k1ar = new arma::mat[2];
  yfin = new arma::mat[2];
  ytmp = new arma::mat[2];
  k2ar = new arma::mat[2];
  k3ar = new arma::mat[2];
  k4ar = new arma::mat[2];
  
  arma::mat y1(dim,1);
  arma::mat y2(dim,1);
  arma::mat gamma1(dim,dim);
  
  L = 0;
  double  nsteps = 50;
  
  for(int i=0; i<N; i++) {
    //Solving differential equations on intervals:
    double t1 = dd(i,1); 
    double t2 = dd(i,2);
    
    int jj=0;
    for(int ii=3; ii<M; ii+=2) {
      y1(jj,0) = dd(i,ii);
      y2(jj,0) = dd(i,ii+1);
      jj += 1;
    }
    
    double tdiff = t2-t1;
    if(tdiff > 2) {
      nsteps = 2*tdiff;
    }
    
    double h = tdiff/nsteps;
    
    //Integration:
    gamma1.zeros(); // set gamma1 to zero matrix
    double s = h/3.00*(-1.00)*mu_quad_linear(t1, y1, gamma1, fH, f1H, mu0H, thetaH, QH, Q1, gomp);
    double t = t1;
    out[0] = y1;
    out[1] = gamma1;
    double ifactor;
    
    for(int j = 0; j < nsteps; j++) {
      //Runge-Kutta method:
      func1_quad_linear(k1ar, t, out, fH, f1H, aH, bH, QH, thetaH);
      yfin[0] = out[0] + h/6.00*k1ar[0];
      yfin[1] = out[1] + h/6.00*k1ar[1];
      ytmp[0] = out[0] + h/2.00*k1ar[0];
      ytmp[1] = out[1] + h/2.00*k1ar[1];
      
      func1_quad_linear(k2ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      yfin[0] = yfin[0] + h/3.00*k2ar[0];
      yfin[1] = yfin[1] + h/3.00*k2ar[1];
      ytmp[0] = out[0] + h/2.00*k2ar[0];
      ytmp[1] = out[1] + h/2.00*k2ar[1];
      
      func1_quad_linear(k3ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      yfin[0] = yfin[0] + h/3.00*k3ar[0];
      yfin[1] = yfin[1] + h/3.00*k3ar[1];
      ytmp[0] = out[0] + h*k3ar[0];
      ytmp[1] = out[1] + h*k3ar[1];
      
      func1_quad_linear(k4ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      out[0] = yfin[0] + h/6.00*k4ar[0];
      out[1] = yfin[1] + h/6.00*k4ar[1];
      
      t = t + h;
      
      //Integration:
      if (j == nsteps-1) {
        ifactor = 1.00;
      } else {
        if (((j % 2) == 0) && (j != 0)) {
          ifactor = 2.00;
        } else {
          ifactor = 4.00;
        }
      }
      
      s = s + ifactor*h/3.00*(-1.00)*mu_quad_linear(t,out[0],out[1], fH, f1H, mu0H, thetaH, QH, Q1, gomp);
    }
    
    arma::mat m2 = out[0];
    arma::mat gamma2 = out[1];
    double pi = 3.141592654;
    
    if(dd(i,0) == 0) { 
      //arma::mat exp = -0.50*dim*log(2.00*pi*det(gamma2)) - 0.50*(m2-y2).t()*pinv(gamma2,ptol)*(m2-y2);
      arma::mat exp = log(pow(2.00*pi, -0.5*dim)*pow(det(gamma2), -0.5)) - 0.50*(m2-y2).t()*pinv(gamma2,ptol)*(m2-y2);
      //arma::mat exp = -0.50*dim*log(2.00*pi*det(gamma2)) - 0.50*(m2-y2).t()*inv(gamma2)*(m2-y2); // inv() fails very ofter
      L += s + exp(0,0);
    } else {
      //double logprobi = log(1.00 - exp(-1.00*mu(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH)));
      double logprobi = log(mu_quad_linear(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH, Q1, gomp));
      L += s + logprobi;
    }
  }
  
  delete[] out;
  delete[] k1ar;
  delete[] yfin;
  delete[] ytmp;
  delete[] k2ar;
  delete[] k3ar;
  delete[] k4ar;
  
  return(Rcpp::wrap(L));
}

RcppExport SEXP complikMD_linear(SEXP dat, 
                                           SEXP n, 
                                           SEXP m, 
                                           SEXP ah, 
                                           SEXP f1h, 
                                           SEXP qh, 
                                           SEXP bh, 
                                           SEXP fh, 
                                           SEXP mu0h, 
                                           SEXP thetah, 
                                           SEXP k, 
                                           SEXP pinv_tol, 
                                           SEXP gomp_) {
  
  long N = as<long>(n); // Number of rows
  long M = as<long>(m); // Number of columns
  
  arma::mat aH = as<arma::mat>(ah);
  arma::mat f1H = as<arma::mat>(f1h);
  arma::mat QH = as<arma::mat>(qh);
  arma::mat bH = as<arma::mat>(bh);
  arma::mat fH = as<arma::mat>(fh);
  double mu0H = as<double>(mu0h); 
  double thetaH  = as<double>(thetah);
  int dim = as<int>(k);
  bool gomp = as<bool>(gomp_);
  
  //Actual data set
  arma::mat dd = as<arma::mat>(dat);  
  double ptol = as<double>(pinv_tol);
  double L; // Likelihood
  
  //End of data loading
  arma::mat *out, *k1ar, *yfin, *ytmp, *k2ar, *k3ar, *k4ar;
  out = new arma::mat[2];
  k1ar = new arma::mat[2];
  yfin = new arma::mat[2];
  ytmp = new arma::mat[2];
  k2ar = new arma::mat[2];
  k3ar = new arma::mat[2];
  k4ar = new arma::mat[2];
  
  arma::mat y1(dim,1);
  arma::mat y2(dim,1);
  arma::mat gamma1(dim,dim);
  
  L = 0;
  double  nsteps = 50;
  
  for(int i=0; i<N; i++) {
    //Solving differential equations on intervals:
    double t1 = dd(i,1); 
    double t2 = dd(i,2);
    
    int jj=0;
    for(int ii=3; ii<M; ii+=2) {
      y1(jj,0) = dd(i,ii);
      y2(jj,0) = dd(i,ii+1);
      jj += 1;
    }
    
    double tdiff = t2-t1;
    if(tdiff > 2) {
      nsteps = 2*tdiff;
    }
    
    double h = tdiff/nsteps;
    
    //Integration:
    gamma1.zeros(); // set gamma1 to zero matrix
    double s = h/3.00*(-1.00)*mu_linear(t1, y1, gamma1, fH, f1H, mu0H, thetaH, QH, gomp);
    double t = t1;
    out[0] = y1;
    out[1] = gamma1;
    double ifactor;
    
    for(int j = 0; j < nsteps; j++) {
      //Runge-Kutta method:
      func1_quad_linear(k1ar, t, out, fH, f1H, aH, bH, QH, thetaH);
      yfin[0] = out[0] + h/6.00*k1ar[0];
      yfin[1] = out[1] + h/6.00*k1ar[1];
      ytmp[0] = out[0] + h/2.00*k1ar[0];
      ytmp[1] = out[1] + h/2.00*k1ar[1];
      
      func1_quad_linear(k2ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      yfin[0] = yfin[0] + h/3.00*k2ar[0];
      yfin[1] = yfin[1] + h/3.00*k2ar[1];
      ytmp[0] = out[0] + h/2.00*k2ar[0];
      ytmp[1] = out[1] + h/2.00*k2ar[1];
      
      func1_quad_linear(k3ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      yfin[0] = yfin[0] + h/3.00*k3ar[0];
      yfin[1] = yfin[1] + h/3.00*k3ar[1];
      ytmp[0] = out[0] + h*k3ar[0];
      ytmp[1] = out[1] + h*k3ar[1];
      
      func1_quad_linear(k4ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      out[0] = yfin[0] + h/6.00*k4ar[0];
      out[1] = yfin[1] + h/6.00*k4ar[1];
      
      t = t + h;
      
      //Integration:
      if (j == nsteps-1) {
        ifactor = 1.00;
      } else {
        if (((j % 2) == 0) && (j != 0)) {
          ifactor = 2.00;
        } else {
          ifactor = 4.00;
        }
      }
      
      s = s + ifactor*h/3.00*(-1.00)*mu_linear(t,out[0],out[1], fH, f1H, mu0H, thetaH, QH, gomp);
    }
    
    arma::mat m2 = out[0];
    arma::mat gamma2 = out[1];
    double pi = 3.141592654;
    
    if(dd(i,0) == 0) { 
      //arma::mat exp = -0.50*dim*log(2.00*pi*det(gamma2)) - 0.50*(m2-y2).t()*pinv(gamma2,ptol)*(m2-y2);
      arma::mat exp = log(pow(2.00*pi, -0.5*dim)*pow(det(gamma2), -0.5)) - 0.50*(m2-y2).t()*pinv(gamma2,ptol)*(m2-y2);
      //arma::mat exp = -0.50*dim*log(2.00*pi*det(gamma2)) - 0.50*(m2-y2).t()*inv(gamma2)*(m2-y2); // inv() fails very ofter
      L += s + exp(0,0);
    } else {
      //double logprobi = log(1.00 - exp(-1.00*mu(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH)));
      double logprobi = log(mu_linear(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH, gomp));
      L += s + logprobi;
    }
  }
  
  delete[] out;
  delete[] k1ar;
  delete[] yfin;
  delete[] ytmp;
  delete[] k2ar;
  delete[] k3ar;
  delete[] k4ar;
  
  return(Rcpp::wrap(L));
}




