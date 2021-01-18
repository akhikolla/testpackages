

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

/*============FUNCTION DEFINITIONS*/
bool gomp = false;

long double mu(double t, arma::mat y1, arma::mat gamma1, arma::mat fH, arma::mat f1H, double mu0H, double thetaH, arma::mat QH) {
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
  mu = mu0Ht + (hfH.t()*QH)*hfH + arma::sum((QH_gamma1).diag());
  
  return (long double)mu(0,0);
}



//Calculating m (y[1]) & gamma(y[2]):
void func1(arma::mat *res, double t, arma::mat *y, arma::mat fH, arma::mat f1H, arma::mat aH, 
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


//Calculating m (y[1]) & gamma(y[2]):
std::vector<arma::mat> func1(double t, std::vector<arma::mat> y, arma::mat fH, arma::mat f1H, arma::mat aH, 
           arma::mat bH, arma::mat QH, double theta) {
  
  arma::mat hfH, hf1H, dy1, dy2;
  //hfH = fH.t() - y[0]; 
  //hf1H = f1H.t() - y[0];
  hfH = fH - y[0]; 
  hf1H = f1H - y[0];
  std::vector<arma::mat> res;
  res.resize(2);
  res[0] = -1.00 * (aH*hf1H) + 2.00 * ((y[1]*Q(t, QH, theta))*hfH);
  res[1] = aH*y[1] + y[1]*aH.t() + bH*bH.t() - 2.00 * ((y[1]*Q(t, QH, theta))*y[1]);
  //std::cout << "theta: " << theta << " QH:" << QH << " y[0]:" << y[0] << " y[1]:" << y[1] << " Q:" << Q(t, QH, theta) << " res0:" << res[0] << " res1:" << res[1]<< std::endl;
  return res;
}



arma::mat Q(double t, arma::mat QH, double theta) {
  arma::mat Q;
  /*if(gomp) {
    Q = QH*exp(theta*t);
  } else {
    Q = QH;
  }*/
  Q = QH;
  return Q;
}


RcppExport SEXP complikMD(SEXP dat, SEXP n, SEXP m, SEXP ah, SEXP f1h, SEXP qh, SEXP bh, SEXP fh, SEXP mu0h, SEXP thetah, SEXP k, SEXP pinv_tol, SEXP gomp_, SEXP logmu0_) {
    
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
    gomp = as<bool>(gomp_);
    bool logmu0 = as<bool>(logmu0_);
    
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
      double s = h/3.00*(-1.00)*mu(t1, y1, gamma1, fH, f1H, logmu0 ? log(mu0H): mu0H, thetaH, QH);
      double t = t1;
      out[0] = y1;
      out[1] = gamma1;
      double ifactor;
      
      for(int j = 0; j < nsteps; j++) {
         //Runge-Kutta method:
         func1(k1ar, t, out, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = out[0] + h/6.00*k1ar[0];
         yfin[1] = out[1] + h/6.00*k1ar[1];
         ytmp[0] = out[0] + h/2.00*k1ar[0];
         ytmp[1] = out[1] + h/2.00*k1ar[1];
         
         func1(k2ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = yfin[0] + h/3.00*k2ar[0];
         yfin[1] = yfin[1] + h/3.00*k2ar[1];
         ytmp[0] = out[0] + h/2.00*k2ar[0];
         ytmp[1] = out[1] + h/2.00*k2ar[1];
         
         func1(k3ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
         yfin[0] = yfin[0] + h/3.00*k3ar[0];
         yfin[1] = yfin[1] + h/3.00*k3ar[1];
         ytmp[0] = out[0] + h*k3ar[0];
         ytmp[1] = out[1] + h*k3ar[1];
         
         func1(k4ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
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
        
        s = s + ifactor*h/3.00*(-1.00)*mu(t,out[0],out[1], fH, f1H, logmu0 ? log(mu0H): mu0H, thetaH, QH);
        
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
        double logprobi = mu(t2, m2, gamma2, fH, f1H, logmu0 ? log(mu0H): mu0H, thetaH, QH);
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

// Simulation routine
RcppExport SEXP simCont(SEXP n, SEXP ah, SEXP f1h, SEXP qh, SEXP fh, SEXP bh, SEXP mu0h, 
                        SEXP thetah, SEXP tstart_, SEXP ystart_, SEXP tend_, SEXP k_, 
                        SEXP dt_, SEXP sd_, SEXP nobs_, SEXP gomp_) {
  
  long N = as<long>(n); // Number of individuals
  
  arma::mat aH = as<arma::mat>(ah);
  arma::mat f1H = as<arma::mat>(f1h);
  arma::mat QH = as<arma::mat>(qh);
  arma::mat bH = as<arma::mat>(bh);
  arma::mat fH = as<arma::mat>(fh);
  double mu0H = as<double>(mu0h); 
  double thetaH  = as<double>(thetah);
  int dim = as<int>(k_);
  gomp = as<bool>(gomp_);
  
  Rcpp::NumericVector sd = Rcpp::NumericVector(sd_);
  int nobs = 0;
  if(!Rf_isNull(nobs_)) {
    nobs = as<int>(nobs_);
  }
  
  
  // Supporting variables
  
  arma::mat *out, *k1ar, *yfin, *ytmp, *k2ar, *k3ar, *k4ar;
  out = new arma::mat[2];
  k1ar = new arma::mat[2];
  yfin = new arma::mat[2];
  ytmp = new arma::mat[2];
  k2ar = new arma::mat[2];
  k3ar = new arma::mat[2];
  k4ar = new arma::mat[2];
  
  /*
  std::vector<arma::mat> out; out.resize(2);
  std::vector<arma::mat> k1ar; k1ar.resize(2);
  std::vector<arma::mat> yfin; yfin.resize(2);
  std::vector<arma::mat> ytmp; ytmp.resize(2);
  std::vector<arma::mat> k2ar; k2ar.resize(2);
  std::vector<arma::mat> k3ar; k3ar.resize(2);
  std::vector<arma::mat> k4ar; k4ar.resize(2);
  */
  
  arma::mat y1(dim,1);
  arma::mat y2(dim,1);
  arma::mat gamma1(dim,dim);
  
  double  nsteps = 50;
  
  double tstart  = as<double>(tstart_);
  arma::mat ystart = as<arma::mat>(ystart_);
  
  double tend  = as<double>(tend_);
  //End of data loading
  double t1;
  double t2;
  double dt=as<double>(dt_);
  bool new_person = false;
  
  std::vector< std::vector<double> > data;
  double S;
  int n_observ;
  
  for(int i=0; i<N; i++) {
    // Starting point
    //t1 = Rcpp::runif(1, tstart, tend)[0];
    t1 = Rcpp::runif(1, tstart, tstart+10)[0];
    
    t2 = t1 + dt + Rcpp::runif(1, 0.0, 1)[0];
    
    for(int ii=0; ii < dim; ii++) {
      y1(ii,0) = Rcpp::rnorm(1, ystart(ii,0), sd[ii])[0];
    }
    
    new_person = false;
    n_observ = 0;
    while(new_person == false) {
      
      double tdiff = t2-t1;
      
      double h = tdiff/nsteps;
      
      gamma1.zeros(); // set gamma1 to zero matrix
      double s;
      s = h/3.00*(-1.00)*mu(t1, y1, gamma1, fH, f1H, mu0H, thetaH, QH);
      
      double t = t1;
      out[0] = y1;
      out[1] = gamma1;
      double ifactor;
      
      for(int j = 0; j < nsteps; j++) {
        //Runge-Kutta method:
        func1(k1ar, t, out, fH, f1H, aH, bH, QH, thetaH);
        yfin[0] = out[0] + h/6.00*k1ar[0];
        yfin[1] = out[1] + h/6.00*k1ar[1];
        ytmp[0] = out[0] + h/2.00*k1ar[0];
        ytmp[1] = out[1] + h/2.00*k1ar[1];
        
        func1(k2ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
        yfin[0] = yfin[0] + h/3.00*k2ar[0];
        yfin[1] = yfin[1] + h/3.00*k2ar[1];
        ytmp[0] = out[0] + h/2.00*k2ar[0];
        ytmp[1] = out[1] + h/2.00*k2ar[1];
        
        func1(k3ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
        yfin[0] = yfin[0] + h/3.00*k3ar[0];
        yfin[1] = yfin[1] + h/3.00*k3ar[1];
        ytmp[0] = out[0] + h*k3ar[0];
        ytmp[1] = out[1] + h*k3ar[1];
        
        func1(k4ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
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
        
        s = s + ifactor*h/3.00*(-1.00)*mu(t,out[0],out[1], fH, f1H, mu0H, thetaH, QH);
        
      }
      
      arma::mat m2 = out[0];
      arma::mat gamma2 = out[1];
      
      S = exp(s);
      
      double xi = 0; // case (0 - alive, 1 - dead) indicator
      
      if(S > Rcpp::runif(1, 0.0, 1.0)[0]) {
        new_person = false;
        
        xi = 0; // case (0 - alive, 1 - dead) indicator
        
        // New y2:
        for(int ii = 0; ii < dim; ii++) {
          y2(ii,0) = Rcpp::rnorm(1, m2(ii,0), sqrt(gamma2(ii,ii)))[0];
        }
        
      } 
      else {
        xi = 1;
        y2 = arma::mat(dim,1);
        
        for(int ii=0; ii<dim; ii++) {
          y2(ii,0) = NumericVector::get_na();
        }
        
        new_person = true;
      }
      
      
      std::vector<double> row; row.resize(4+2*dim);
      row[0] = i; row[1] = xi; row[2] = t1; row[3] = t2;
      
      int jj=0;
      for(int ii=4; ii<(4 + 2*dim-1); ii+=2) {
        row[ii] = y1(jj,0);
        row[ii+1] = y2(jj,0);
        jj += 1;
      }
      data.push_back(row);
      n_observ += 1;
      if(n_observ == nobs) {
        new_person = true;
      }
      
      if(new_person == false) {
        y1 = y2;
        t1 = t2;
        t2 = t1 + dt + Rcpp::runif(1, 0.0, 1)[0];
        if(t2 > tend) {
          new_person = true;
          break;
        }
      }
    }
  }
  
  delete[] out;
  delete[] k1ar;
  delete[] yfin;
  delete[] ytmp;
  delete[] k2ar;
  delete[] k3ar;
  delete[] k4ar;
  
  arma::mat res(data.size(), 4+2*dim);
  for(int i=0; i<data.size(); i++) {
    for(int j=0; j<4+2*dim; j++) {
      res(i, j) = data[i][j];
    }
  }
  
  return(Rcpp::wrap(res));
}

// Likelihood for GenSPM (see Arbeev et al., "Genetic model for longitudinal studies of aging, health 
// and longevity and its potential application incomplete data" 2011)
RcppExport SEXP complik_gen(SEXP dat, SEXP n, SEXP m, 
                            SEXP ah, SEXP al, 
                            SEXP f1h, SEXP f1l, 
                            SEXP qh, SEXP ql, 
                            SEXP bh, SEXP bl, 
                            SEXP fh, SEXP fl, 
                            SEXP mu0h, SEXP mu0l, 
                            SEXP thetah, SEXP thetal, 
                            SEXP p_, 
                            SEXP nc_, SEXP nnc_, 
                            SEXP k, 
                            SEXP pinv_tol,
                            SEXP gomp_) {
  
  long N = as<long>(n); // Number of rows
  long M = as<long>(m); // Number of columns
  
  arma::mat aH = as<arma::mat>(ah); 
  arma::mat aL = as<arma::mat>(al);
  arma::mat f1H = as<arma::mat>(f1h); 
  arma::mat f1L = as<arma::mat>(f1l);
  arma::mat QH = as<arma::mat>(qh); 
  arma::mat QL = as<arma::mat>(ql); 
  arma::mat bH = as<arma::mat>(bh); 
  arma::mat bL = as<arma::mat>(bl);
  arma::mat fH = as<arma::mat>(fh); 
  arma::mat fL = as<arma::mat>(fl); 
  double mu0H = as<double>(mu0h); double mu0L = as<double>(mu0l); 
  double thetaH  = as<double>(thetah); double thetaL  = as<double>(thetal);
  int dim = as<int>(k);
  int nc = as<int>(nc_);
  int nnc = as<int>(nnc_);
  gomp = as<bool>(gomp_);
  
  //Actual data set
  arma::mat dd = as<arma::mat>(dat);  
  double ptol = as<double>(pinv_tol);
  double L; // Likelihood
  double p = as<double>(p_);
  
  //End of data loading
  /*arma::mat *out, *k1ar, *yfin, *ytmp, *k2ar, *k3ar, *k4ar;
  out = new arma::mat[2];
  k1ar = new arma::mat[2]; k1ar->resize(2,1);
  yfin = new arma::mat[2]; yfin->resize(2,1);
  ytmp = new arma::mat[2]; ytmp->resize(2,1);
  k2ar = new arma::mat[2]; k2ar->resize(2,1);
  k3ar = new arma::mat[2]; k3ar->resize(2,1);
  k4ar = new arma::mat[2]; k4ar->resize(2,1);
  */
  std::vector<arma::mat> out; out.resize(2);
  std::vector<arma::mat> k1ar; k1ar.resize(2);
  std::vector<arma::mat> yfin; yfin.resize(2);
  std::vector<arma::mat> ytmp; ytmp.resize(2);
  std::vector<arma::mat> k2ar; k2ar.resize(2);
  std::vector<arma::mat> k3ar; k3ar.resize(2);
  std::vector<arma::mat> k4ar; k4ar.resize(2);
  
  arma::mat y1(dim,1);
  arma::mat y2(dim,1);
  arma::mat gamma1(dim,dim);
  
  L = 0.0;
  double  nsteps = 20;
  
  int G = 0; // Genetic variable
  
  for(int i=0; i<N; i++) {
    G = dd(i, 3);
    //Solving differential equations on intervals:
    double t1 = dd(i,1); 
    double t2 = dd(i,2);
    
    int jj=0;
    for(int ii=4; ii<M; ii+=2) { //for(int ii=3; ii<M; ii+=2) {
      y1(jj,0) = dd(i,ii);
      y2(jj,0) = dd(i,ii+1);
      jj += 1;
    }
    
    double tdiff = t2-t1;
    /*if(tdiff > 2) {
      nsteps = 2*tdiff;
    }*/
    
    double h = tdiff/nsteps;
    
    //Integration:
    gamma1.zeros(); // set gamma1 to zero matrix
    double s;
    if(G == 1) {
      s = h/3.00*(-1.00)*mu(t1, y1, gamma1, fH, f1H, mu0H, thetaH, QH);
    } else {
      s = h/3.00*(-1.00)*mu(t1, y1, gamma1, fL, f1L, mu0L, thetaL, QL);
    }
    
    double t = t1;
    out[0] = y1;
    out[1] = gamma1;
    double ifactor;
    
    for(int j = 0; j < nsteps; j++) {
      //Runge-Kutta method:
      if(G == 1) {
        k1ar = func1(t, out, fH, f1H, aH, bH, QH, thetaH);
      } else {
        k1ar = func1(t, out, fL, f1L, aL, bL, QL, thetaL);
      }
      yfin[0] = out[0] + h/6.00*k1ar[0];
      yfin[1] = out[1] + h/6.00*k1ar[1];
      ytmp[0] = out[0] + h/2.00*k1ar[0];
      ytmp[1] = out[1] + h/2.00*k1ar[1];
      
      if(G == 1) {
        k2ar = func1(t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      } else {
        k2ar = func1(t, ytmp, fL, f1L, aL, bL, QL, thetaL);
      }
      yfin[0] = yfin[0] + h/3.00*k2ar[0];
      yfin[1] = yfin[1] + h/3.00*k2ar[1];
      ytmp[0] = out[0] + h/2.00*k2ar[0];
      ytmp[1] = out[1] + h/2.00*k2ar[1];
      
      if(G == 1) {
        k3ar = func1(t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      } else {
        k3ar = func1(t, ytmp, fL, f1L, aL, bL, QL, thetaL);
      }
      yfin[0] = yfin[0] + h/3.00*k3ar[0];
      yfin[1] = yfin[1] + h/3.00*k3ar[1];
      ytmp[0] = out[0] + h*k3ar[0];
      ytmp[1] = out[1] + h*k3ar[1];
      
      if(G == 1) {
        k4ar = func1(t, ytmp, fH, f1H, aH, bH, QH, thetaH);
      } else {
        k4ar = func1(t, ytmp, fL, f1L, aL, bL, QL, thetaL);
      }
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
      
      if(G == 1) {
        s = s + ifactor*h/3.00*(-1.00)*mu(t,out[0],out[1], fH, f1H, mu0H, thetaH, QH);
      } else {
        s = s + ifactor*h/3.00*(-1.00)*mu(t,out[0],out[1], fL, f1L, mu0L, thetaL, QL);
      }
      
    }
    
    arma::mat m2 = out[0];
    arma::mat gamma2 = out[1];
    double pi = 3.141592654;
    
    if(dd(i,0) == 0) { 
      long double det_gamma = det(gamma2);
      if(det_gamma < 0) {
        det_gamma = 1e-16;
      }
      //arma::mat exp = -0.50*dim*log(2.00*pi*pow(det(gamma2), -0.5)) - 0.50*(m2-y2).t()*pinv(gamma2,ptol)*(m2-y2);
      arma::mat exp = log(pow(2.00*pi, -0.5*dim)*pow(det_gamma, -0.5)) - 0.50*(m2-y2).t()*pinv(gamma2,ptol)*(m2-y2);
      //arma::mat exp = -0.50*dim*log(2.00*pi*det(gamma2)) - 0.50*(m2-y2).t()*inv(gamma2)*(m2-y2); // inv() fails very ofter
      
      L += s + exp(0,0); 
      
    } else {
      double logprobi; 
      if(G == 1) {
        //logprobi = log(1.00 - exp(-1.00*mu(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH)));
        logprobi = log(mu(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH));
        
      } else {
        //logprobi = log(1.00 - exp(-1.00*mu(t2, m2, gamma2, fL, f1L, mu0L, thetaL, QL)));
        logprobi = log(mu(t2, m2, gamma2, fL, f1L, mu0L, thetaL, QL));
       
      }
      
      L += s + logprobi;
      
    }
    
  }
  
  /*delete[] out;
  delete[] k1ar;
  delete[] yfin;
  delete[] ytmp;
  delete[] k2ar;
  delete[] k3ar;
  delete[] k4ar;
  */
  //L = pow(p, nc)*pow(1-p, nnc)*L;
  L = log(pow(p, nc)*pow(1-p, nnc)) + L;
  
  return(Rcpp::wrap(L));
}

// Simulation routine for genetic data
RcppExport SEXP simGenCont(SEXP n, 
                           SEXP ah, SEXP al, 
                           SEXP f1h, SEXP f1l, 
                           SEXP qh, SEXP ql, 
                           SEXP fh, SEXP fl,
                           SEXP bh, SEXP bl,
                           SEXP mu0h, SEXP mu0l,
                           SEXP thetah, SEXP thetal,
                           SEXP p0_,
                           SEXP tstart_, SEXP ystart_, SEXP tend_, 
                           SEXP k_, SEXP dt_, SEXP sd_, 
                           SEXP _genmode, SEXP gomp_, SEXP nobs_) {
  
  long N = as<long>(n); // Number of individuals
  
  arma::mat aH = as<arma::mat>(ah); 
  arma::mat aL = as<arma::mat>(al);
  arma::mat f1H = as<arma::mat>(f1h); 
  arma::mat f1L = as<arma::mat>(f1l);
  arma::mat QH = as<arma::mat>(qh); 
  arma::mat QL = as<arma::mat>(ql); 
  arma::mat bH = as<arma::mat>(bh); 
  arma::mat bL = as<arma::mat>(bl);
  arma::mat fH = as<arma::mat>(fh); 
  arma::mat fL = as<arma::mat>(fl); 
  
  long double mu0H = as<double>(mu0h); 
  double mu0L = as<double>(mu0l); 
  double thetaH  = as<double>(thetah); 
  double thetaL  = as<double>(thetal);
  int dim = as<int>(k_);
  Rcpp::NumericVector sd = Rcpp::NumericVector(sd_);
  double p0 = as<double>(p0_);
  int genmode = as<int>(_genmode);
  
  gomp = as<bool>(gomp_);
  int nobs = as<int>(nobs_);
  
  // Supporting variables
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
  
  double  nsteps = 20;
  
  Rcpp::NumericVector tstart = Rcpp::NumericVector(tstart_);
  arma::mat ystart = as<arma::mat>(ystart_);
  
  double tend  = as<double>(tend_);
  //End of data loading
  double t1;
  double t2;
  long double dt=as<long double>(dt_);
  bool new_person = false;
  
  std::vector< std::vector<double> > data;
  double S;
  int n_observ;
  
  for(int i=0; i<N; i++) {
    // Starting point
    //t1 = Rcpp::runif(1, tstart, tend)[0];
    //t1 = R::runif(tstart, tstart+10);
    //t1 = Rcpp::runif(1, tstart, tstart+10)[0];
    
    if(tstart.size() == 1) {
      t1 = Rcpp::runif(1, tstart[0], tend)[0]; //Starting time
    } else if(tstart.size() == 2) {
      t1 = Rcpp::runif(1,tstart[0], tstart[1])[0]; //Starting time
    }
    
    //t2 = t1 + dt + Rcpp::runif(1, 0.0,1)[0]; 
    t2 = t1 + Rcpp::runif(1,-1.0*dt/10,dt/10)[0] + dt;
    
    new_person = false;
    
    for(int ii=0; ii < dim; ii++) {
      y1(ii,0) = Rcpp::rnorm(1, ystart(ii,0), sd[ii])[0];
    }
    
    new_person = false;
    
    int G = 0;
    if(Rcpp::runif(1, 0.0, 1.00)[0] < p0) {
      G = 1;
    }
    //G = Rcpp::rbinom(1,1,p0)[0];
    
    
    n_observ = 0;
    
    while(new_person == false) {
      
      double tdiff = t2-t1;
      
      double h = tdiff/nsteps;
      
      gamma1.zeros(); // set gamma1 to zero matrix
      double s;
      if(G == 1) {
        s = h/3.00*(-1.00)*mu(t1, y1, gamma1, fH, f1H, mu0H, thetaH, QH);
      } else {
        s = h/3.00*(-1.00)*mu(t1, y1, gamma1, fL, f1L, mu0L, thetaL, QL);
      }
      
      double t = t1;
      out[0] = y1;
      out[1] = gamma1;
      double ifactor;
      
      for(int j = 0; j < nsteps; j++) {
        //Runge-Kutta method:
        if(G == 1) {
          func1(k1ar, t, out, fH, f1H, aH, bH, QH, thetaH);
        } else {
          func1(k1ar, t, out, fL, f1L, aL, bL, QL, thetaL);
        }
        
        yfin[0] = out[0] + h/6.00*k1ar[0];
        yfin[1] = out[1] + h/6.00*k1ar[1];
        ytmp[0] = out[0] + h/2.00*k1ar[0];
        ytmp[1] = out[1] + h/2.00*k1ar[1];
        
        if(G == 1) {
          func1(k2ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
        } else {
          func1(k2ar, t, ytmp, fL, f1L, aL, bL, QL, thetaL);
        }
        
        yfin[0] = yfin[0] + h/3.00*k2ar[0];
        yfin[1] = yfin[1] + h/3.00*k2ar[1];
        ytmp[0] = out[0] + h/2.00*k2ar[0];
        ytmp[1] = out[1] + h/2.00*k2ar[1];
        
        if(G == 1) {
          func1(k3ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
        } else {
          func1(k3ar, t, ytmp, fL, f1L, aL, bL, QL, thetaL);
        }
        
        yfin[0] = yfin[0] + h/3.00*k3ar[0];
        yfin[1] = yfin[1] + h/3.00*k3ar[1];
        ytmp[0] = out[0] + h*k3ar[0];
        ytmp[1] = out[1] + h*k3ar[1];
        
        if(G == 1) {
          func1(k4ar, t, ytmp, fH, f1H, aH, bH, QH, thetaH);
        } else {
          func1(k4ar, t, ytmp, fL, f1L, aL, bL, QL, thetaL);
        }
        
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
        
        if(G == 1) {
          s = s + ifactor*h/3.00*(-1.00)*mu(t,out[0],out[1], fH, f1H, mu0H, thetaH, QH);
        } else {
          s = s + ifactor*h/3.00*(-1.00)*mu(t,out[0],out[1], fL, f1L, mu0L, thetaL, QL);
        }
      }
      
      arma::mat m2 = out[0];
      arma::mat gamma2 = out[1];
      
      S = exp(s);
      
      double xi = 0; // case (0 - alive, 1 - dead) indicator
      
      if(S > Rcpp::runif(1, 0.0, 1.0)[0]) { //if(S > R::runif(0.0, 1.0)) {
        new_person = false;
        xi = 0; // case (0 - alive, 1 - dead) indicator
        // New y2:
        for(int ii = 0; ii < dim; ii++) {
          y2(ii,0) = Rcpp::rnorm(1, m2(ii,0), sqrt(gamma2(ii,ii)))[0];
        }
      } 
      else {
        xi = 1;
        y2 = arma::mat(dim,1);
        
        for(int ii=0; ii<dim; ii++) {
          y2(ii,0) = NumericVector::get_na();
        }
        
        new_person = true;
      }
      
      std::vector<double> row; 
      int jj=0;
      if(genmode == 1) {
        row.resize(5+2*dim);
        row[0] = i; row[1] = xi; row[2] = t1; row[3] = t2;
        row[4] = G;
        for(int ii=5; ii<(5 + 2*dim); ii+=2) {
          row[ii] = y1(jj,0);
          row[ii+1] = y2(jj,0);
          jj += 1;
        }
      } else {
        row.resize(4+2*dim);
        row[0] = i; row[1] = xi; row[2] = t1; row[3] = t2;
        for(int ii=4; ii<(4 + 2*dim); ii+=2) {
          row[ii] = y1(jj,0);
          row[ii+1] = y2(jj,0);
          jj += 1;
        }
      }
      
      data.push_back(row);
      
      n_observ += 1;
      if(n_observ == nobs) {
        new_person = true;
      }
      
      if(new_person == false) {
        y1 = y2;
        t1 = t2;
        //t2 = t1 + dt + Rcpp::runif(1, 0.0,1)[0]; 
        t2 = t1 + Rcpp::runif(1,-1.0*dt/10,dt/10)[0] + dt;
        if(t2 > tend) {
          new_person = true;
          break;
        }
      }
    }
  }
  
  delete[] out;
  delete[] k1ar;
  delete[] yfin;
  delete[] ytmp;
  delete[] k2ar;
  delete[] k3ar;
  delete[] k4ar;
  
  if(genmode == 1) {
    arma::mat res(data.size(), 5+2*dim);
    for(int i=0; i<data.size(); i++) {
      for(int j=0; j<5+2*dim; j++) {
        res(i, j) = data[i][j];
      }
    }
    return(Rcpp::wrap(res));
  } else {
    arma::mat res(data.size(), 4+2*dim);
    for(int i=0; i<data.size(); i++) {
      for(int j=0; j<4+2*dim; j++) {
        res(i, j) = data[i][j];
      }
    }
    return(Rcpp::wrap(res));
  }
}

// Likelihood for GenSPM (see Arbeev et al., "Genetic model for longitudinal studies of aging, health 
// and longevity and its potential application incomplete data" 2011)
// This particular function estimates likelihood from a non-genetic group
RcppExport SEXP complikGenNonGenetic(SEXP dat, SEXP n, SEXP m, 
                            SEXP ah, SEXP al, 
                            SEXP f1h, SEXP f1l, 
                            SEXP qh, SEXP ql, 
                            SEXP bh, SEXP bl, 
                            SEXP fh, SEXP fl, 
                            SEXP mu0h, SEXP mu0l, 
                            SEXP thetah, SEXP thetal, 
                            SEXP p_, 
                            SEXP k, 
                            SEXP pinv_tol, SEXP gomp_) {
  
  long N = as<long>(n); // Number of rows
  long M = as<long>(m); // Number of columns
  
  // Carriers: H, non-carriers: L
  arma::mat aH = as<arma::mat>(ah); 
  arma::mat aL = as<arma::mat>(al);
  arma::mat f1H = as<arma::mat>(f1h); 
  arma::mat f1L = as<arma::mat>(f1l);
  arma::mat QH = as<arma::mat>(qh); 
  arma::mat QL = as<arma::mat>(ql); 
  arma::mat bH = as<arma::mat>(bh); 
  arma::mat bL = as<arma::mat>(bl);
  arma::mat fH = as<arma::mat>(fh); 
  arma::mat fL = as<arma::mat>(fl); 
  gomp = as<bool>(gomp_);
  long double mu0H = as<long double>(mu0h); 
  long double mu0L = as<long double>(mu0l); 
  long double thetaH  = as<long double>(thetah); 
  long double thetaL  = as<long double>(thetal);
  int dim = as<int>(k);
  
  //Actual data set
  arma::mat dd = as<arma::mat>(dat);  
  long double ptol = as<long double>(pinv_tol);
  long double L, L1, L0; // Likelihoods for carriers (1) and non-carriers (0)
  long double p = as<long double>(p_);
  //===============================End of data loading========================================//
  
  // Carriers
  /*
  arma::mat *out1, *k1ar1, *yfin1, *ytmp1, *k2ar1, *k3ar1, *k4ar1;
  out1 = new arma::mat[2];
  k1ar1 = new arma::mat[2];
  yfin1 = new arma::mat[2];
  ytmp1 = new arma::mat[2];
  k2ar1 = new arma::mat[2];
  k3ar1 = new arma::mat[2];
  k4ar1 = new arma::mat[2];
  */
  std::vector<arma::mat> out1; out1.resize(2);
  std::vector<arma::mat> k1ar1; k1ar1.resize(2);
  std::vector<arma::mat> yfin1; yfin1.resize(2);
  std::vector<arma::mat> ytmp1; ytmp1.resize(2);
  std::vector<arma::mat> k2ar1; k2ar1.resize(2);
  std::vector<arma::mat> k3ar1; k3ar1.resize(2);
  std::vector<arma::mat> k4ar1; k4ar1.resize(2);
  //=====================//
  // Non-carriers
  /*
  arma::mat *out0, *k1ar0, *yfin0, *ytmp0, *k2ar0, *k3ar0, *k4ar0;
  out0 = new arma::mat[2];
  k1ar0 = new arma::mat[2];
  yfin0 = new arma::mat[2];
  ytmp0 = new arma::mat[2];
  k2ar0 = new arma::mat[2];
  k3ar0 = new arma::mat[2];
  k4ar0 = new arma::mat[2];
  */
  std::vector<arma::mat> out0; out0.resize(2);
  std::vector<arma::mat> k1ar0; k1ar0.resize(2);
  std::vector<arma::mat> yfin0; yfin0.resize(2);
  std::vector<arma::mat> ytmp0; ytmp0.resize(2);
  std::vector<arma::mat> k2ar0; k2ar0.resize(2);
  std::vector<arma::mat> k3ar0; k3ar0.resize(2);
  std::vector<arma::mat> k4ar0; k4ar0.resize(2);
  
  //=====================//
  
  arma::mat y1(dim,1);
  arma::mat y2(dim,1);
  arma::mat gamma1(dim,dim);
  
  L = 0.0;
  long double  nsteps = 20;
  
  for(int i=0; i<N; i++) {
    //Solving differential equations on intervals:
    long double t1 = dd(i,1); 
    long double t2 = dd(i,2);
    
    int jj=0;
    for(int ii=3; ii<M; ii+=2) { // ii=0: id, ii=1: t1, ii=2: t2, ii=3: y1, ii=4: y2 (y.next)
      y1(jj,0) = dd(i,ii);
      y2(jj,0) = dd(i,ii+1);
      jj += 1;
    }
    
    long double tdiff = t2-t1;
    /*if(tdiff > 2) {
      nsteps = 2*tdiff;
    }*/
    
    long double h = tdiff/nsteps;
    
    //Integration:
    gamma1.zeros(); // set gamma1 to zero matrix
    long double s1, s0; // For carriers (1) and non-carriers (0)
    
    // Carriers
    s1 = h/3.00*(-1.00)*mu(t1, y1, gamma1, fH, f1H, mu0H, thetaH, QH);
    gamma1.zeros();
    s0 = h/3.00*(-1.00)*mu(t1, y1, gamma1, fL, f1L, mu0L, thetaL, QL);
    
    long double t = t1;
    
    // Carriers
    out1[0] = y1;
    out1[1] = gamma1;
    // Non-carriers
    out0[0] = y1;
    out0[1] = gamma1;
    
    long double ifactor;
    long double L1Y, L1Q, L0Y, L0Q;
    L1Y = L1Q = L0Y = L0Q = 1.0;
    for(int j = 0; j < nsteps; j++) {
      //Runge-Kutta method:
      // Carriers
      //func1(k1ar1, t, out1, fH, f1H, aH, bH, QH, thetaH);
      k1ar1 = func1(t, out1, fH, f1H, aH, bH, QH, thetaH);
      // Non-carriers
      //func1(k1ar0, t, out0, fL, f1L, aL, bL, QL, thetaL);
      k1ar0 = func1(t, out0, fL, f1L, aL, bL, QL, thetaL);
      // Carriers
      yfin1[0] = out1[0] + h/6.00*k1ar1[0];
      yfin1[1] = out1[1] + h/6.00*k1ar1[1];
      ytmp1[0] = out1[0] + h/2.00*k1ar1[0];
      ytmp1[1] = out1[1] + h/2.00*k1ar1[1];
      // Non-carriers
      yfin0[0] = out0[0] + h/6.00*k1ar0[0];
      yfin0[1] = out0[1] + h/6.00*k1ar0[1];
      ytmp0[0] = out0[0] + h/2.00*k1ar0[0];
      ytmp0[1] = out0[1] + h/2.00*k1ar0[1];
      
      // Carriers
      //func1(k2ar1, t, ytmp1, fH, f1H, aH, bH, QH, thetaH);
      k2ar1 = func1(t, ytmp1, fH, f1H, aH, bH, QH, thetaH);
      // Non-carriers
      //func1(k2ar0, t, ytmp0, fL, f1L, aL, bL, QL, thetaL);
      k2ar0 = func1(t, ytmp0, fL, f1L, aL, bL, QL, thetaL);
      // Carriers
      yfin1[0] = yfin1[0] + h/3.00*k2ar1[0];
      yfin1[1] = yfin1[1] + h/3.00*k2ar1[1];
      ytmp1[0] = out1[0] + h/2.00*k2ar1[0];
      ytmp1[1] = out1[1] + h/2.00*k2ar1[1];
      // Non-carriers
      yfin0[0] = yfin0[0] + h/3.00*k2ar0[0];
      yfin0[1] = yfin0[1] + h/3.00*k2ar0[1];
      ytmp0[0] = out0[0] + h/2.00*k2ar0[0];
      ytmp0[1] = out0[1] + h/2.00*k2ar0[1];
      
      // Carriers
      //func1(k3ar1, t, ytmp1, fH, f1H, aH, bH, QH, thetaH);
      k3ar1 = func1(t, ytmp1, fH, f1H, aH, bH, QH, thetaH);
      // Non-carriers
      //func1(k3ar0, t, ytmp0, fL, f1L, aL, bL, QL, thetaL);
      k3ar0 = func1(t, ytmp0, fL, f1L, aL, bL, QL, thetaL);
      // Carriers
      yfin1[0] = yfin1[0] + h/3.00*k3ar1[0];
      yfin1[1] = yfin1[1] + h/3.00*k3ar1[1];
      ytmp1[0] = out1[0] + h*k3ar1[0];
      ytmp1[1] = out1[1] + h*k3ar1[1];
      // Non-carriers
      yfin0[0] = yfin0[0] + h/3.00*k3ar0[0];
      yfin0[1] = yfin0[1] + h/3.00*k3ar0[1];
      ytmp0[0] = out0[0] + h*k3ar0[0];
      ytmp0[1] = out0[1] + h*k3ar0[1];
      
      // Carriers
      //func1(k4ar1, t, ytmp1, fH, f1H, aH, bH, QH, thetaH);
      k4ar1 = func1(t, ytmp1, fH, f1H, aH, bH, QH, thetaH);
      // Non-carriers
      //func1(k4ar0, t, ytmp0, fL, f1L, aL, bL, QL, thetaL);
      k4ar0 = func1(t, ytmp0, fL, f1L, aL, bL, QL, thetaL);
      // Carriers
      out1[0] = yfin1[0] + h/6.00*k4ar1[0];
      out1[1] = yfin1[1] + h/6.00*k4ar1[1];
      // Non-carriers
      out0[0] = yfin0[0] + h/6.00*k4ar0[0];
      out0[1] = yfin0[1] + h/6.00*k4ar0[1];
      
      //=================================================//
      
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
      
      // Carriers
      s1 = s1 + ifactor*h/3.00*(-1.00)*mu(t,out1[0],out1[1], fH, f1H, mu0H, thetaH, QH);
      // Non-carriers
      s0 = s0 + ifactor*h/3.00*(-1.00)*mu(t,out0[0],out0[1], fL, f1L, mu0L, thetaL, QL);
      
      
    }
    
    // Carriers
    arma::mat m21 = out1[0];
    arma::mat gamma21 = out1[1];
    // Non-carriers
    arma::mat m20 = out0[0];
    arma::mat gamma20 = out0[1];
    
    long double pi = 3.141592654;
    
    bool firstRow = false;
    bool lastRow = false;
    
    if(i == 0) {
     firstRow = true; 
    } else if(dd(i,1) != dd(i-1,2)) {
      firstRow = true;
    }
    
    if((i < N-1) && (dd(i,2) != dd(i+1,1))){
      lastRow = true;
    }else if(i == N-1) {
      lastRow = true;
    }
    
    arma::mat exp1;
    arma::mat exp0;
    
    
    if(firstRow) {
      L1 = 1.0; L0 = 1.0;
      if(dd(i,0) == 0) {
        exp1 = exp(-0.5*(m21-y2).t()*pinv(gamma21,ptol)*(m21-y2));
        exp0 = exp(-0.5*(m20-y2).t()*pinv(gamma20,ptol)*(m20-y2));
        
        L1 = L1*pow(2.00*pi,-0.5*dim)*1/sqrt(det(gamma21))*exp1(0,0)*exp(s1);
        L0 = L0*pow(2.00*pi,-0.5*dim)*1/sqrt(det(gamma20))*exp0(0,0)*exp(s0);
      } else {
        L1 = L1*mu(t2, m21, gamma21, fH, f1H, mu0H, thetaH, QH)*exp(s1);
        L0 = L0*mu(t2, m20, gamma20, fL, f1L, mu0L, thetaL, QL)*exp(s0);
      }
    } else {
      if(dd(i,0) == 0) {
        exp1 = exp(-0.5*(m21-y2).t()*pinv(gamma21,ptol)*(m21-y2));
        exp0 = exp(-0.5*(m20-y2).t()*pinv(gamma20,ptol)*(m20-y2));
        
        L1 = L1*pow(2.00*pi,-0.5*dim)*1/sqrt(det(gamma21))*exp1(0,0)*exp(s1);
        L0 = L0*pow(2.00*pi,-0.5*dim)*1/sqrt(det(gamma20))*exp0(0,0)*exp(s0);
      } else {
        L1 = L1*mu(t2, m21, gamma21, fH, f1H, mu0H, thetaH, QH)*exp(s1);
        L0 = L0*mu(t2, m20, gamma20, fL, f1L, mu0L, thetaL, QL)*exp(s0);
      }
    }
    
    if(lastRow) {
      
      L += log(p*L1 + (1-p)*L0);
    }
     
     
  }
  
  // Cleaning up
  /*delete[] out1;
  delete[] k1ar1;
  delete[] yfin1;
  delete[] ytmp1;
  delete[] k2ar1;
  delete[] k3ar1;
  delete[] k4ar1;
  //
  delete[] out0;
  delete[] k1ar0;
  delete[] yfin0;
  delete[] ytmp0;
  delete[] k2ar0;
  delete[] k3ar0;
  delete[] k4ar0;
  */
  return(Rcpp::wrap(L));
}

