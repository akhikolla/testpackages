#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// Deconvolution Kernel grid evaluations
RcppExport SEXP Ku_sec_order(SEXP ngrid_, SEXP delta_, SEXP h1_, SEXP sigU_, SEXP t_, SEXP dt_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  int ngrid = as<int>(ngrid_);
  double delta = as<double>(delta_);
  NumericVector h1(h1_);
  double sigU = as<double>(sigU_);
  NumericVector t(t_);
  double dt = as<double>(dt_);
  int nh1 = h1.size();
  
  // results to save
  NumericMatrix Ku0(ngrid, nh1);
  NumericMatrix Ku1(ngrid, nh1);
  NumericMatrix Ku2(ngrid, nh1);
  NumericVector FKt = FK_sec_order(t);
  NumericVector FKt1 = FK1_sec_order(t);
  NumericVector FKt2 = FK2_sec_order(t);
  // deconvolution kernel;
  for(int k=0; k<nh1; ++k){
    NumericVector FUtinv = FuLapinv(t/h1[k], sigU);
    NumericVector FKt_FUt = FKt*FUtinv;
    NumericVector FKt1_FUt = FKt1*FUtinv;
    NumericVector FKt2_FUt = FKt2*FUtinv;
    for(int i=0; i<ngrid; ++i){
      Ku0(i,k) = Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt_FUt)*dt/(2.0*PI);
      //if(Ku0(i,k)<0) { Ku0(i,k)=0; break; }
      Ku1(i,k) = -Rcpp::sum(Rcpp::sin((i+0.0)*delta*t)*FKt1_FUt)*dt/(2.0*PI);
      Ku2(i,k) = -Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt2_FUt)*dt/(2.0*PI);
    }
  }
  return List::create(Named("Ku0")=Ku0,
                      Named("Ku1")=Ku1,
                      Named("Ku2")=Ku2);
  END_RCPP
}
RcppExport SEXP Ku0_sec_order(SEXP ngrid_, SEXP delta_, SEXP h1_, SEXP sigU_, SEXP t_, SEXP dt_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  int ngrid = as<int>(ngrid_);
  double delta = as<double>(delta_);
  NumericVector h1(h1_);
  double sigU = as<double>(sigU_);
  NumericVector t(t_);
  double dt = as<double>(dt_);
  int nh1 = h1.size();
  
  // results to save
  NumericMatrix Ku0(ngrid, nh1);
  NumericVector FKt = FK_sec_order(t);
  // deconvolution kernel;
  for(int k=0; k<nh1; ++k){
    NumericVector FKt_FUt = FKt*FuLapinv(t/h1[k], sigU);
    for(int i=0; i<ngrid; ++i){
      Ku0(i,k) = Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt_FUt)*dt/(2.0*PI);
      //if(Ku0(i,k)<0) { Ku0(i,k)=0; break; }
    }
  }
  return List::create(Named("Ku0")=Ku0);
  END_RCPP
}

// estimate the mode curve when error is laplace, local constant
RcppExport SEXP LCfitModeRegLap2( SEXP x_, SEXP y_, SEXP yindx_, SEXP W_, SEXP Y_, SEXP dt_, SEXP t_,
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
  
  // deconvolution kernel
  double delta=0.005;
  double max1 = std::abs(Rcpp::max(W)-Rcpp::min(x));
  double max2 = std::abs(Rcpp::max(W)-Rcpp::max(x));
  double max3 = std::abs(Rcpp::min(W)-Rcpp::min(x));
  double max4 = std::abs(Rcpp::min(W)-Rcpp::max(x));
  int ngrid = round(std::max(std::max(max1, max2), std::max(max3, max4))/h1/delta)+1;
  //Rprintf( "iscan = %d\n", ngrid );
  NumericVector Ku0(ngrid);
  for(int i=0; i<ngrid; ++i){
    Ku0[i] = Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt_FUt)*dt/(2.0*PI);
    //if(Ku0[i]<0) { Ku0[i]=0; break; }
  }
  
  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      int indx = round(std::abs(W[i]-x[j])/h1/delta);
      // if(indx>ngrid) indx=ngrid;
      Ku0ij(i,j) = Ku0[indx];
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
          //std::abs(KG_tot)<1e-10
          //KG_tot<=1e-10
          newy=NA_REAL;
          break;
        }else{
          newy = YKG_tot/KG_tot;
          err_now = std::abs(newy-oldy);
          iter_now++;
        }
      }
      //if(iter_now==max_iterations) newy=NA_REAL;
      if((iter_now==max_iterations)&&(err_now > (eps*10))) newy=NA_REAL;
      ym[jj]=newy;
    }
  }
  return List::create(Named("mode")=ym);
  END_RCPP
}

// estimate the mode curve when error is laplace, local linear estimator
RcppExport SEXP LLfitModeRegLap2( SEXP x_, SEXP y_, SEXP yindx_, SEXP W_, SEXP Y_, SEXP dt_, SEXP t_,
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
  
  // deconvolution kernel
  double delta=0.005;
  double max1 = std::abs(Rcpp::max(W)-Rcpp::min(x));
  double max2 = std::abs(Rcpp::max(W)-Rcpp::max(x));
  double max3 = std::abs(Rcpp::min(W)-Rcpp::min(x));
  double max4 = std::abs(Rcpp::min(W)-Rcpp::max(x));
  int ngrid = round(std::max(std::max(max1, max2), std::max(max3, max4))/h1/delta)+1;
  // Rprintf( "iscan = %d\n", ngrid );
  NumericVector Ku0(ngrid);
  NumericVector Ku1(ngrid);
  NumericVector Ku2(ngrid);
  for(int i=0; i<ngrid; ++i){
    Ku0[i] = Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt_FUt)*dt/(2.0*PI);
    //if(Ku0[i]<0) { Ku0[i]=0; break; }
    Ku1[i] = -Rcpp::sum(Rcpp::sin((i+0.0)*delta*t)*FKt1_FUt)*dt/(2.0*PI);
    Ku2[i] = -Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt2_FUt)*dt/(2.0*PI);
  }
  
  for(int j=0; j<nx; ++j){
    S1=0; S2=0;
    for(int i=0; i<n; ++i){
      int indx = round(std::abs(W[i]-x[j])/h1/delta);
      //if(indx>ngrid) indx=ngrid;
      K0i[i] = Ku0[indx];
      if((W[i]-x[j])<0) {
        K1i[i] = 0.0-Ku1[indx];
      }else{
        K1i[i] = Ku1[indx];
      }
      S1 += K1i[i];
      S2 += Ku2[indx];
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
      //if(iter_now==max_iterations) newy=NA_REAL;
      if((iter_now==max_iterations)&&(err_now > (eps*10))) newy=NA_REAL;
      ym[jj]=newy;
    }
  }
  return List::create(Named("mode")=ym);
  END_RCPP
}
