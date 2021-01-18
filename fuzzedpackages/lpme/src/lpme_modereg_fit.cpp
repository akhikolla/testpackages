#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// estimate the mode curve without measurement error, local constant
RcppExport SEXP LCfitModeReg(SEXP x_, SEXP y_, SEXP yindx_, SEXP X_, SEXP Y_, 
                              SEXP h1_, SEXP h2_, SEXP max_iterations_, SEXP eps_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector x(x_);
  NumericVector y(y_);
  IntegerVector yindx(yindx_);
  NumericVector X(X_);
  NumericVector Y(Y_);
  double h1 = as<double>(h1_);
  double h2 = as<double>(h2_);
  int max_iterations = as<int>(max_iterations_);
  double eps = as<double>(eps_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  NumericVector ym(ny);
  
  // temp variables
  NumericMatrix Ku0ij(n,nx);
  double KGj = 0;
  double KG_tot=0;
  double YKG_tot=0;
  double oldy=0;
  double newy=0;
  int iter_now;
  double err_now;
  
  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      Ku0ij(i,j) = exp( -0.5*std::pow(((X[i]-x[j])/h1), 2) );
    }
  }
  
  for(int j=0; j<nx; ++j){
    //Rprintf( "iscan = %d\n", j+1 );
    R_CheckUserInterrupt();
    int ind1 = yindx[j];
    int ind2 = yindx[j+1]-1;
    for(int jj=ind1; jj<=ind2; ++jj){
      newy = y[jj];
      err_now=1e10;
      iter_now=0;
      while((iter_now < max_iterations)&&(err_now > eps)){
        //R_CheckUserInterrupt();
        YKG_tot=0;
        KG_tot=0;
        oldy = newy;
        for(int i=0; i<n; ++i){
          KGj = Ku0ij(i,j)*exp( -0.5*std::pow(((newy-Y[i])/h2), 2) );
          KG_tot += KGj;
          YKG_tot += Y[i]*KGj;
        }
        if((KG_tot)<1e-10){
          newy=NA_REAL;
          break;
        }else{
          newy = YKG_tot/KG_tot;
          err_now = std::abs(newy-oldy);
          iter_now++;
        }
      }
      if((iter_now==max_iterations)&&(err_now > (eps*10))) newy=NA_REAL;
      ym[jj]=newy;
    }
  }
  return List::create(Named("mode")=ym);
  END_RCPP
}

// estimate the mode curve, local linear estimator
RcppExport SEXP LLfitModeReg( SEXP x_, SEXP y_, SEXP yindx_, SEXP X_, SEXP Y_, 
                               SEXP h1_, SEXP h2_, SEXP max_iterations_, SEXP eps_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector x(x_);
  NumericVector y(y_);
  IntegerVector yindx(yindx_);
  NumericVector X(X_);
  NumericVector Y(Y_);
  double h1 = as<double>(h1_);
  double h2 = as<double>(h2_);
  int max_iterations = as<int>(max_iterations_);
  double eps = as<double>(eps_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  NumericVector ym(ny);
  
  // temp variables
  NumericMatrix Ku0ij(n,nx);
  NumericVector K0i(n);
  NumericVector K1i(n);
  double S1=0, S2=0;
  double KGj = 0;
  double KG_tot=0;
  double YKG_tot=0;
  double oldy=0;
  double newy=0;
  int iter_now;
  double err_now;
  
  for(int j=0; j<nx; ++j){
    S1=0; S2=0;
    for(int i=0; i<n; ++i){
      K0i[i] = exp( -0.5*std::pow(((X[i]-x[j])/h1), 2) );
      K1i[i] = ((X[i]-x[j])/h1)*K0i[i];
      S1 += K1i[i];
      S2 += ((X[i]-x[j])/h1)*K1i[i];
    }
    for(int i=0; i<n; ++i){
      Ku0ij(i,j) = (K0i[i]*S2 - K1i[i]*S1);
    }
  }
  
  for(int j=0; j<nx; ++j){
    //Rprintf( "iscan = %d\n", j+1 );
    R_CheckUserInterrupt();
    int ind1 = yindx[j];
    int ind2 = yindx[j+1]-1;
    for(int jj=ind1; jj<=ind2; ++jj){
      newy = y[jj];
      err_now=1e10;
      iter_now=0;
      while((iter_now < max_iterations)&&(err_now > eps)){
        //R_CheckUserInterrupt();
        YKG_tot=0;
        KG_tot=0;
        oldy = newy;
        for(int i=0; i<n; ++i){
          KGj = Ku0ij(i,j)*exp( -0.5*std::pow(((newy-Y[i])/h2), 2) );
          KG_tot += KGj;
          YKG_tot += Y[i]*KGj;
        }
        if((KG_tot)<1e-10){
          newy=NA_REAL;
          break;
        }else{
          newy = YKG_tot/KG_tot;
          err_now = std::abs(newy-oldy);
          iter_now++;
        }
      }
      if((iter_now==max_iterations)&&(err_now > (eps*10))) newy=NA_REAL;
      ym[jj]=newy;
    }
  }
  return List::create(Named("mode")=ym);
  END_RCPP
}

// estimate the mode curve when error is laplace, local constant
RcppExport SEXP LCfitModeRegLap( SEXP x_, SEXP y_, SEXP yindx_, SEXP W_, SEXP Y_, SEXP dt_, SEXP t_,
                                  SEXP sigU_, SEXP h1_, SEXP h2_, SEXP max_iterations_, SEXP eps_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector x(x_);
  NumericVector y(y_);
  IntegerVector yindx(yindx_);
  NumericVector W(W_);
  NumericVector Y(Y_);
  double dt = as<double>(dt_);
  NumericVector t(t_);
  double sigU = as<double>(sigU_);
  double h1 = as<double>(h1_);
  double h2 = as<double>(h2_);
  int max_iterations = as<int>(max_iterations_);
  double eps = as<double>(eps_);
  int nx = x.size();
  int ny = y.size();
  int n = W.size();
  NumericVector FKt_FUt = FK_sec_order(t)*FuLapinv(t/h1, sigU);
  
  // results to save 
  NumericVector ym(ny);
  
  // temp variables
  NumericMatrix Ku0ij(n,nx);
  double KGj = 0;
  double KG_tot=0;
  double YKG_tot=0;
  double oldy=0;
  double newy=0;
  int iter_now;
  double err_now;
  
  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      Ku0ij(i,j) = Rcpp::sum(Rcpp::cos((W[i]-x[j])/h1*t)*FKt_FUt)*dt;
    }
  }
  
  for(int j=0; j<nx; ++j){
    //Rprintf( "iscan = %d\n", j+1 );
    R_CheckUserInterrupt();
    int ind1 = yindx[j];
    int ind2 = yindx[j+1]-1;
    for(int jj=ind1; jj<=ind2; ++jj){
      newy = y[jj];
      err_now=1e10;
      iter_now=0;
      while((iter_now < max_iterations)&&(err_now > eps)){
        //R_CheckUserInterrupt();
        YKG_tot=0;
        KG_tot=0;
        oldy = newy;
        for(int i=0; i<n; ++i){
          KGj = Ku0ij(i,j)*exp( -0.5*std::pow(((newy-Y[i])/h2), 2) );
          KG_tot += KGj;
          YKG_tot += Y[i]*KGj;
        }
        if((KG_tot)<1e-10){
          newy=NA_REAL;
          break;
        }else{
          newy = YKG_tot/KG_tot;
          err_now = std::abs(newy-oldy);
          iter_now++;
        }
      }
      if((iter_now==max_iterations)&&(err_now > (eps*10))) newy=NA_REAL;
      ym[jj]=newy;
    }
  }
  return List::create(Named("mode")=ym);
  END_RCPP
}

// estimate the mode curve when error is laplace, local linear estimator
RcppExport SEXP LLfitModeRegLap( SEXP x_, SEXP y_, SEXP yindx_, SEXP W_, SEXP Y_, SEXP dt_, SEXP t_,
                                  SEXP sigU_, SEXP h1_, SEXP h2_, SEXP max_iterations_, SEXP eps_) {
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector x(x_);
  NumericVector y(y_);
  IntegerVector yindx(yindx_);
  NumericVector W(W_);
  NumericVector Y(Y_);
  double dt = as<double>(dt_);
  NumericVector t(t_);
  double sigU = as<double>(sigU_);
  double h1 = as<double>(h1_);
  double h2 = as<double>(h2_);
  int max_iterations = as<int>(max_iterations_);
  double eps = as<double>(eps_);
  int nx = x.size();
  int ny = y.size();
  int n = W.size();
  NumericVector FUtinv = FuLapinv(t/h1, sigU);
  NumericVector FKt_FUt = FK_sec_order(t)*FUtinv;
  NumericVector FKt1_FUt = FK1_sec_order(t)*FUtinv;
  NumericVector FKt2_FUt = FK2_sec_order(t)*FUtinv;
  
  
  // results to save 
  NumericVector ym(ny);
  
  // temp variables
  NumericMatrix Ku0ij(n,nx);
  NumericVector K0i(n);
  NumericVector K1i(n);
  double S1=0, S2=0;
  double KGj = 0;
  double KG_tot=0;
  double YKG_tot=0;
  double oldy=0;
  double newy=0;
  int iter_now;
  double err_now;
  
  for(int j=0; j<nx; ++j){
    S1=0; S2=0;
    for(int i=0; i<n; ++i){
      K0i[i] = Rcpp::sum(Rcpp::cos((W[i]-x[j])/h1*t)*FKt_FUt)*dt;
      K1i[i] = -Rcpp::sum(Rcpp::sin((W[i]-x[j])/h1*t)*FKt1_FUt)*dt;
      S1 += K1i[i];
      S2 += (-Rcpp::sum(Rcpp::cos((W[i]-x[j])/h1*t)*FKt2_FUt)*dt);
    }
    for(int i=0; i<n; ++i){
      Ku0ij(i,j) = (K0i[i]*S2 - K1i[i]*S1);
    }
  }
  
  for(int j=0; j<nx; ++j){
    //Rprintf( "iscan = %d\n", j+1 );
    R_CheckUserInterrupt();
    int ind1 = yindx[j];
    int ind2 = yindx[j+1]-1;
    for(int jj=ind1; jj<=ind2; ++jj){
      newy = y[jj];
      err_now=1e10;
      iter_now=0;
      while((iter_now < max_iterations)&&(err_now > eps)){
        //R_CheckUserInterrupt();
        YKG_tot=0;
        KG_tot=0;
        oldy = newy;
        for(int i=0; i<n; ++i){
          KGj = Ku0ij(i,j)*exp( -0.5*std::pow(((newy-Y[i])/h2), 2) );
          KG_tot += KGj;
          YKG_tot += Y[i]*KGj;
        }
        if((KG_tot)<1e-10){
          newy=NA_REAL;
          break;
        }else{
          newy = YKG_tot/KG_tot;
          err_now = std::abs(newy-oldy);
          iter_now++;
        }
      }
      if((iter_now==max_iterations)&&(err_now > (eps*10))) newy=NA_REAL;
      ym[jj]=newy;
    }
  }
  return List::create(Named("mode")=ym);
  END_RCPP
}
