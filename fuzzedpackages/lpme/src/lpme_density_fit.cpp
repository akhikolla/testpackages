#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// Univariate kernel density without measurement error, Gaussian kernel
RcppExport SEXP fitDensityGauK(SEXP X_, SEXP x_, SEXP h1_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector x(x_);
  double h1=as<double>(h1_);
  int nx = x.size();
  int n = X.size();
  
  // results to save 
  NumericVector res(nx);
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    double Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x +=exp( -0.5*std::pow(((X[i]-x[ii])/h1), 2) )/sqrt2pi;
    }
    res[ii] = Ku0x/nh1;
  }
  return List::create(Named("fit")=res);
  END_RCPP
}

// Univariate kernel density without measurement error, second order kernel
RcppExport SEXP fitDensitySecK(SEXP X_, SEXP x_, SEXP h1_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector x(x_);
  double h1=as<double>(h1_);
  int nx = x.size();
  int n = X.size();
  
  // results to save 
  NumericVector res(nx);
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    double Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x +=K_sec_order((X[i]-x[ii])/h1);
    }
    res[ii] = Ku0x/nh1;
  }
  return List::create(Named("fit")=res);
  END_RCPP
}

// Univariate kernel density with Laplace measurement error, second order kernel
RcppExport SEXP fitDensitySecKLap(SEXP W_, SEXP x_, SEXP h1_, SEXP dt_, SEXP t_, SEXP sigU_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector x(x_);
  double h1=as<double>(h1_);
  int nx = x.size();
  int n = W.size();
  double dt = as<double>(dt_);
  NumericVector t(t_);
  double sigU = as<double>(sigU_);
  NumericVector FKt_FUt = FK_sec_order(t)*FuLapinv(t/h1, sigU);
  
  // results to save 
  NumericVector res(nx);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  
  // deconvolution kernel
  double delta=0.01;
  double max1 = std::abs(Rcpp::max(W)-Rcpp::min(x));
  double max2 = std::abs(Rcpp::max(W)-Rcpp::max(x));
  double max3 = std::abs(Rcpp::min(W)-Rcpp::min(x));
  double max4 = std::abs(Rcpp::min(W)-Rcpp::max(x));
  int ngrid = round(std::max(std::max(max1, max2), std::max(max3, max4))/h1/delta)+1;
  if(ngrid<(n*nx)){
    NumericVector Ku0(ngrid);
    for(int i=0; i<ngrid; ++i){
      Ku0[i] = Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt_FUt)*dt/(2.0*PI);
    }
    for(int i=0; i<n; ++i){
      for(int j=0; j<nx; ++j){
        int indx = round(std::abs(W[i]-x[j])/h1/delta);
        Ku0ij(i,j) = Ku0[indx];
      }
    }
  }else{
    for(int i=0; i<n; ++i){
      for(int j=0; j<nx; ++j){
        Ku0ij(i,j) = Rcpp::sum(Rcpp::cos((W[i]-x[j])/h1*t)*FKt_FUt)*dt/(2.0*PI);
      }
    }
  }
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    double Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x +=Ku0ij(i,ii);
    }
    if(Ku0x<0){
      res[ii] = 0.0;
    }else{
      res[ii] = Ku0x/nh1;
    }
  }
  return List::create(Named("fit")=res);
  END_RCPP
}

// Multivariate kernel density without measurement error, Gaussian kernel
RcppExport SEXP fitDensityGauK2(SEXP X_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  NumericMatrix res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  
  // kernel for X
  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      Ku0ij(i,j) = exp( -0.5*std::pow(((X[i]-x[j])/h1), 2) )/sqrt2pi;
    }
  }
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    for(int jj=0; jj<ny; ++jj){
      double Ku0xGy=0;
      for(int i=0; i<n; ++i){
        Ku0xGy += Ku0ij(i, ii)*exp( -0.5*std::pow(((Y[i]-y[jj])/h2), 2) )/sqrt2pi;
      }
      res(ii,jj) = Ku0xGy/nh1/h2;
    }
  }
  return List::create(Named("fit")=res);
  END_RCPP
}

// Multivariate kernel density without measurement error, second order kernel
RcppExport SEXP fitDensitySecK2(SEXP X_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  NumericMatrix res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  
  // kernel for X
  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      Ku0ij(i,j) = K_sec_order((X[i]-x[j])/h1);
    }
  }
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    for(int jj=0; jj<ny; ++jj){
      double Ku0xGy=0;
      for(int i=0; i<n; ++i){
        Ku0xGy += Ku0ij(i, ii)*exp( -0.5*std::pow(((Y[i]-y[jj])/h2), 2) )/sqrt2pi;
      }
      res(ii,jj) = Ku0xGy/nh1/h2;
    }
  }
  return List::create(Named("fit")=res);
  END_RCPP
}

// Multivariate kernel density with Laplace measurement error, second order kernel
RcppExport SEXP fitDensitySecKLap2(SEXP W_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_, 
                                   SEXP dt_, SEXP t_, SEXP sigU_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  double dt = as<double>(dt_);
  NumericVector t(t_);
  double sigU = as<double>(sigU_);
  int nx = x.size();
  int ny = y.size();
  int n = W.size();
  NumericVector FKt_FUt = FK_sec_order(t)*FuLapinv(t/h1, sigU);
  
  // results to save 
  NumericMatrix res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  
  // deconvolution kernel
  double delta=0.01;
  double max1 = std::abs(Rcpp::max(W)-Rcpp::min(x));
  double max2 = std::abs(Rcpp::max(W)-Rcpp::max(x));
  double max3 = std::abs(Rcpp::min(W)-Rcpp::min(x));
  double max4 = std::abs(Rcpp::min(W)-Rcpp::max(x));
  int ngrid = round(std::max(std::max(max1, max2), std::max(max3, max4))/h1/delta)+1;
  if(ngrid<(n*nx)){
    NumericVector Ku0(ngrid);
    for(int i=0; i<ngrid; ++i){
      Ku0[i] = Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt_FUt)*dt/(2.0*PI);
    }
    for(int i=0; i<n; ++i){
      for(int j=0; j<nx; ++j){
        int indx = round(std::abs(W[i]-x[j])/h1/delta);
        Ku0ij(i,j) = Ku0[indx];
      }
    }
  }else{
    for(int i=0; i<n; ++i){
      for(int j=0; j<nx; ++j){
        Ku0ij(i,j) = Rcpp::sum(Rcpp::cos((W[i]-x[j])/h1*t)*FKt_FUt)*dt/(2.0*PI);
      }
    }
  }
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    for(int jj=0; jj<ny; ++jj){
      double Ku0xGy=0;
      for(int i=0; i<n; ++i){
        Ku0xGy += Ku0ij(i, ii)*exp( -0.5*std::pow(((Y[i]-y[jj])/h2), 2) )/sqrt2pi;
      }
      if(Ku0xGy<0){
        res(ii,jj) = 0.0;
      }else{
        res(ii,jj) = Ku0xGy/nh1/h2;
      }
    }
  }
  return List::create(Named("fit")=res);
  END_RCPP
}

