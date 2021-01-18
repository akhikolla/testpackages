#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// second-order Kernel K
RcppExport SEXP Kern_2nd_order(SEXP x_){
  NumericVector x(x_);
  int n = x.size();
  NumericVector res(n);
  for(int i=0; i<n; ++i){
    double xi = std::abs(x[i]);
    if(xi<0.2){
      res[i] = 0.1455068+0.0000996*xi+ -0.0084387*std::pow(xi,2);
    }else{
      res[i] = 48.0*std::cos(xi)/(PI*std::pow(xi,4))*(1.0-15.0/std::pow(xi,2)) - 
        144*std::sin(xi)/(PI*std::pow(xi,5))*(2.0-5.0/std::pow(xi,2));
    }
  }
  return wrap(res);
}

// estimate conditional density without measurement error, local constant
RcppExport SEXP LCfitDensityReg(SEXP X_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_, SEXP mx_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  NumericVector mx(mx_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  arma::mat res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  double Ku0x=0;
  double Ku0xGy=0;
  
  // kernel for X
  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      Ku0ij(i,j) = exp( -0.5*std::pow(((X[i]-x[j])/h1), 2) )/sqrt2pi;
    }
  }
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    R_CheckUserInterrupt();
    Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x += Ku0ij(i,ii);
    }
    if((Ku0x/nh1)<1e-10){
      (res.row(ii)).fill(0);
    }else{
      for(int jj=0; jj<ny; ++jj){
        Ku0xGy=0;
        for(int i=0; i<n; ++i){
          Ku0xGy += Ku0ij(i, ii)*exp( -0.5*std::pow(((Y[i]-y[jj]+mx[ii])/h2), 2) )/sqrt2pi;
        }
        if(Ku0xGy<0){
          res(ii, jj) = 0.0;
        }else{
          res(ii, jj) = Ku0xGy/Ku0x/h2;
        }
      }
    }
  }
  return List::create(Named("fitxy")=res);
  END_RCPP
}

RcppExport SEXP LCfitDensityRegK(SEXP X_, SEXP Y_, SEXP x_, SEXP y_,
                                    SEXP h1_, SEXP h2_, SEXP mx_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  NumericVector mx(mx_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  arma::mat res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  double Ku0x=0;
  double Ku0xGy=0;
  
  // kernel for X
  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      Ku0ij(i,j) = K_sec_order((X[i]-x[j])/h1);
    }
  }
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    R_CheckUserInterrupt();
    Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x += Ku0ij(i,ii);
    }
    if((Ku0x/nh1)<1e-10){
      (res.row(ii)).fill(0);
    }else{
      for(int jj=0; jj<ny; ++jj){
        Ku0xGy=0;
        for(int i=0; i<n; ++i){
          Ku0xGy += Ku0ij(i, ii)*exp( -0.5*std::pow(((Y[i]-y[jj]+mx[ii])/h2), 2) )/sqrt2pi;
        }
        if(Ku0xGy<0){
          res(ii, jj) = 0.0;
        }else{
          res(ii, jj) = Ku0xGy/Ku0x/h2;
        }
      }
    }
  }
  return List::create(Named("fitxy")=res);
  END_RCPP
}

RcppExport SEXP LCfitDensityRegKK(SEXP X_, SEXP Y_, SEXP x_, SEXP y_,
                                    SEXP h1_, SEXP h2_, SEXP mx_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  NumericVector mx(mx_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  arma::mat res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  double Ku0x=0;
  double Ku0xGy=0;
  
  // kernel for X
  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      Ku0ij(i,j) = K_sec_order((X[i]-x[j])/h1);
    }
  }
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    R_CheckUserInterrupt();
    Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x += Ku0ij(i,ii);
    }
    if((Ku0x/nh1)<1e-10){
      (res.row(ii)).fill(0);
    }else{
      for(int jj=0; jj<ny; ++jj){
        Ku0xGy=0;
        for(int i=0; i<n; ++i){
          Ku0xGy += Ku0ij(i, ii)*K_sec_order((Y[i]-y[jj]+mx[ii])/h2);
        }
        if(Ku0xGy<0){
          res(ii, jj) = 0.0;
        }else{
          res(ii, jj) = Ku0xGy/Ku0x/h2;
        }
      }
    }
  }
  return List::create(Named("fitxy")=res);
  END_RCPP
}

// estimate conditional density without measurement error, local linear
RcppExport SEXP LLfitDensityReg(SEXP X_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_, SEXP mx_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  NumericVector mx(mx_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  arma::mat res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  NumericVector K0i(n);
  NumericVector K1i(n);
  double S1=0, S2=0;
  double Ku0x=0;
  double Ku0xGy=0;
  
  // kernel for X
  for(int j=0; j<nx; ++j){
    S1=0; S2=0;
    for(int i=0; i<n; ++i){
      K0i[i] = exp( -0.5*std::pow(((X[i]-x[j])/h1), 2) )/sqrt2pi;
      K1i[i] = ((X[i]-x[j])/h1)*K0i[i];
      S1 += K1i[i];
      S2 += ((X[i]-x[j])/h1)*K1i[i];
    }
    for(int i=0; i<n; ++i){
      Ku0ij(i,j) = (K0i[i]*S2 - K1i[i]*S1);
    }
  }
  
  // start estimate
  double nh1 = std::pow((n+0.0)*h1, 2);
  for(int ii=0; ii<nx; ++ii){
    R_CheckUserInterrupt();
    Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x += Ku0ij(i,ii);
    }
    if((Ku0x/nh1)<1e-10){
      (res.row(ii)).fill(0);
    }else{
      for(int jj=0; jj<ny; ++jj){
        Ku0xGy=0;
        for(int i=0; i<n; ++i){
          Ku0xGy += Ku0ij(i, ii)*exp( -0.5*std::pow(((Y[i]-y[jj]+mx[ii])/h2), 2) )/sqrt2pi;
        }
        if(Ku0xGy<0){
          res(ii, jj) = 0.0;
        }else{
          res(ii, jj) = Ku0xGy/Ku0x/h2;
        }
      }
    }
  }
  return List::create(Named("fitxy")=res);
  END_RCPP
}

// estimate conditional density when error is laplace, local constant, approximated method
RcppExport SEXP LCfitDensityRegLap(SEXP W_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_, SEXP mx_, 
                                   SEXP dt_, SEXP t_, SEXP sigU_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  NumericVector mx(mx_);
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
  arma::mat res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  double Ku0x=0;
  double Ku0xGy=0;
  
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
    R_CheckUserInterrupt();
    Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x += Ku0ij(i,ii);
    }
    if((Ku0x/nh1)<1e-10){
      (res.row(ii)).fill(0);
    }else{
      for(int jj=0; jj<ny; ++jj){
        Ku0xGy=0;
        for(int i=0; i<n; ++i){
          Ku0xGy += Ku0ij(i, ii)*exp( -0.5*std::pow(((Y[i]-y[jj]+mx[ii])/h2), 2) )/sqrt2pi;
        }
        if(Ku0xGy<0){
          res(ii, jj) = 0.0;
        }else{
          res(ii, jj) = Ku0xGy/Ku0x/h2;
        }
      }
    }
  }
  return List::create(Named("fitxy")=res);
  END_RCPP
}

// estimate conditional density when error is laplace, local constant, approximated method
RcppExport SEXP LCfitDensityRegLapLap(SEXP W_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_, SEXP mx_, 
                                   SEXP dt_, SEXP t_, SEXP sigU_, SEXP sigU2_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  NumericVector mx(mx_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  double dt = as<double>(dt_);
  NumericVector t(t_);
  double sigU = as<double>(sigU_);
  double sigU2 = as<double>(sigU2_);
  int nx = x.size();
  int ny = y.size();
  int n = W.size();
  NumericVector FKt_FUt = FK_sec_order(t)*FuLapinv(t/h1, sigU);
  NumericVector FKt_FU2t = FK_sec_order(t)*FuLapinv(t/h2, sigU2);
  
  // results to save 
  arma::mat res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  double Ku0x=0;
  double Ku0xGy=0;
  int kk=0;
  
  // deconvolution kernel
  double delta=0.01;
  double maxW=Rcpp::max(W); double minW=Rcpp::min(W); 
  double maxx=Rcpp::max(x); double minx=Rcpp::min(x); 
  double max1 = std::abs(maxW-minx);
  double max2 = std::abs(maxW-maxx);
  double max3 = std::abs(minW-minx);
  double max4 = std::abs(minW-maxx);
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
  
  // deconvolution kernel for K2
  Rcpp::NumericVector rangeYymx(n*nx*ny);
  kk=0;
  for(int ii=0; ii<nx; ++ii){
    for(int jj=0; jj<ny; ++jj){
      for(int i=0; i<n; ++i){
        rangeYymx[kk]=Y[i]-y[jj]+mx[ii];
        ++kk;
      }
    }
  }
  max1 = std::abs(Rcpp::max(rangeYymx));
  max2 = std::abs(Rcpp::min(rangeYymx));
  int ngrid2 = round(std::max(max1, max2)/h2/delta)+1;
  NumericVector Ku2_0(ngrid2);
  if(ngrid2<(n*nx*ny)){
    for(int i=0; i<ngrid2; ++i){
      Ku2_0[i] = Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt_FU2t)*dt/(2.0*PI);
    }
  }
  
  // start estimate
  double nh1 = (n+0.0)*h1;
  for(int ii=0; ii<nx; ++ii){
    R_CheckUserInterrupt();
    Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x += Ku0ij(i,ii);
    }
    if((Ku0x/nh1)<1e-10){
      (res.row(ii)).fill(0);
    }else{
      for(int jj=0; jj<ny; ++jj){
        Ku0xGy=0;
        for(int i=0; i<n; ++i){
          double Ku2_0Yymx=0;
          if(ngrid2<(n*nx*ny)){
            int indx = round(std::abs(Y[i]-y[jj]+mx[ii])/h2/delta);
            Ku2_0Yymx=Ku2_0[indx];
          }else{
            Ku2_0Yymx=Rcpp::sum(Rcpp::cos((Y[i]-y[jj]+mx[ii])/h2*t)*FKt_FU2t)*dt/(2.0*PI);
          }
          Ku0xGy += Ku0ij(i, ii)*Ku2_0Yymx;
        }
        if(Ku0xGy<0){
          res(ii, jj) = 0.0;
        }else{
          res(ii, jj) = Ku0xGy/Ku0x/h2;
        }
      }
    }
  }
  return List::create(Named("fitxy")=res);
  END_RCPP
}

// estimate conditional density when error is laplace, local linear, approximated method
RcppExport SEXP LLfitDensityRegLap(SEXP W_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_, SEXP mx_, 
                                    SEXP dt_, SEXP t_, SEXP sigU_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  NumericVector mx(mx_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  double dt = as<double>(dt_);
  NumericVector t(t_);
  double sigU = as<double>(sigU_);
  int nx = x.size();
  int ny = y.size();
  int n = W.size();
  NumericVector FUtinv = FuLapinv(t/h1, sigU);
  NumericVector FKt_FUt = FK_sec_order(t)*FUtinv;
  NumericVector FKt1_FUt = FK1_sec_order(t)*FUtinv;
  NumericVector FKt2_FUt = FK2_sec_order(t)*FUtinv;
  
  // results to save 
  arma::mat res(nx, ny);
  
  // temp variable
  NumericMatrix Ku0ij(n,nx);
  NumericVector K0i(n);
  NumericVector K1i(n);
  double S1=0, S2=0;
  double Ku0x=0;
  double Ku0xGy=0;
  
  // deconvolution kernel
  double delta=0.01;
  double max1 = std::abs(Rcpp::max(W)-Rcpp::min(x));
  double max2 = std::abs(Rcpp::max(W)-Rcpp::max(x));
  double max3 = std::abs(Rcpp::min(W)-Rcpp::min(x));
  double max4 = std::abs(Rcpp::min(W)-Rcpp::max(x));
  int ngrid = round(std::max(std::max(max1, max2), std::max(max3, max4))/h1/delta)+1;
  if(ngrid<(n*nx)){
    NumericVector Ku0(ngrid);
    NumericVector Ku1(ngrid);
    NumericVector Ku2(ngrid);
    for(int i=0; i<ngrid; ++i){
      Ku0[i] = Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt_FUt)*dt/(2.0*PI);
      Ku1[i] = -Rcpp::sum(Rcpp::sin((i+0.0)*delta*t)*FKt1_FUt)*dt/(2.0*PI);
      Ku2[i] = -Rcpp::sum(Rcpp::cos((i+0.0)*delta*t)*FKt2_FUt)*dt/(2.0*PI);
    }
    for(int j=0; j<nx; ++j){
      S1=0; S2=0;
      for(int i=0; i<n; ++i){
        int indx = round(std::abs(W[i]-x[j])/h1/delta);
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
  }else{
    for(int j=0; j<nx; ++j){
      S1=0; S2=0;
      for(int i=0; i<n; ++i){
        K0i[i] = Rcpp::sum(Rcpp::cos((W[i]-x[j])/h1*t)*FKt_FUt)*dt/(2.0*PI);
        K1i[i] = -Rcpp::sum(Rcpp::sin((W[i]-x[j])/h1*t)*FKt1_FUt)*dt/(2.0*PI);
        S1 += K1i[i];
        S2 += (-Rcpp::sum(Rcpp::cos((W[i]-x[j])/h1*t)*FKt2_FUt)*dt)/(2.0*PI);
      }
      for(int i=0; i<n; ++i){
        Ku0ij(i,j) = (K0i[i]*S2 - K1i[i]*S1);
      }
    }
  }
  
  // start estimate
  double nh1 = std::pow((n+0.0)*h1, 2);
  for(int ii=0; ii<nx; ++ii){
    R_CheckUserInterrupt();
    Ku0x=0;
    for(int i=0; i<n; ++i){
      Ku0x += Ku0ij(i,ii);
    }
    if((Ku0x/nh1)<1e-10){
      (res.row(ii)).fill(0);
    }else{
      for(int jj=0; jj<ny; ++jj){
        Ku0xGy=0;
        for(int i=0; i<n; ++i){
          Ku0xGy += Ku0ij(i, ii)*exp( -0.5*std::pow(((Y[i]-y[jj]+mx[ii])/h2), 2) )/sqrt2pi;
        }
        if(Ku0xGy<0){
          res(ii, jj) = 0.0;
        }else{
          res(ii, jj) = Ku0xGy/Ku0x/h2;
        }
      }
    }
  }
  return List::create(Named("fitxy")=res);
  END_RCPP
}

// estimate conditional density when error is laplace, local constant, approximated method, HZ method
RcppExport SEXP LCfitDensityRegLap2(SEXP W_, SEXP Y_, SEXP x_, SEXP y_, SEXP h1_, SEXP h2_, SEXP mx_, 
                                    SEXP input_, SEXP output_, SEXP beta_, SEXP beta2_, SEXP mconst_,
                                    SEXP Kinput_, SEXP sigU_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector Y(Y_);
  NumericVector x(x_);
  NumericVector y(y_);
  NumericVector mx(mx_);
  double h1=as<double>(h1_);
  double h2=as<double>(h2_);
  NumericVector input(input_);
  NumericVector output(output_);
  arma::vec mcon = as<arma::vec>(mconst_);
  double beta = as<double>(beta_);
  double beta2 = as<double>(beta2_);
  NumericVector Kinput(Kinput_);
  double sigU = as<double>(sigU_);
  int nx = x.size();
  int ny = y.size();
  int n = W.size();
  
  // results to save 
  arma::mat res(nx, ny);
  
  // temp variable
  int m = input.size();
  int m_mid = m/2 +1; 
  arma::vec fWin(m);
  arma::vec gfWin(m);
  double nh1 = n*h1;
  NumericMatrix Ku0ij(n,m);
  double Ku0x=0;
  double Ku0xGy=0;
  
  // kernel for W
  for(int i=0; i<n; ++i){
    for(int j=0; j<m; ++j){
      int indx = (int)(round((W[i]-input[j])/h1/beta+m_mid));
      Ku0ij(i,j) = ((indx<=m) & (indx>=1))? Kinput[indx-1]:0;
    }
  }
  
  // Find the support for CF of KK in the denominator
  int indexl = (int)(round(-1.0/h1/beta2+m_mid)); 
  int indexu = (int)(round(1.0/h1/beta2+m_mid)); 
  arma::vec support1 = arma::ones<vec>(m);
  for (int i=0; i<(indexl-1); ++i) {support1[i]=0;}
  for (int i=indexu; i<m; ++i) {support1[i]=0;}
  // Find the support for CF of KK in the numerator
  indexl = (int)(round(-2.0/h1/beta2+m_mid)); 
  indexu = (int)(round(2.0/h1/beta2+m_mid)); 
  arma::vec support2 = arma::ones<vec>(m);
  for (int i=0; i<(indexl-1); ++i) {support2[i]=0;}
  for (int i=indexu; i<m; ++i) {support2[i]=0;}
  // FfU
  NumericVector FfU=FuLap(output, sigU);
  arma::vec FfU2(FfU.begin(), FfU.size(), false);
  //start estimating
  for(int jj=0; jj<ny; ++jj){
    R_CheckUserInterrupt();
    for(int ii=0; ii<m; ++ii){
      Ku0x=0;
      for(int i=0; i<n; ++i){
        Ku0x += Ku0ij(i,ii);
      }
      fWin[ii]=Ku0x/nh1;
      Ku0xGy=0;
      for(int i=0; i<n; ++i){
        int indx = (int)(round((Y[i]-y[jj]+mx[ii])/h2/beta+m_mid));
        double Ku0y = ((indx<=m) & (indx>=1))? Kinput[indx-1]:0;
        //double Ku0y = exp( -0.5*std::pow(((Y[i]-y[jj]+mx[ii])/h2), 2) )/sqrt2pi;
        Ku0xGy += Ku0ij(i, ii)*Ku0y;
      }
      gfWin[ii]=Ku0xGy/nh1/h2;
    }
    // FFT for fW
    arma::cx_vec FfW = beta*mcon%arma::fft( mcon%fWin )%support1;
    // FFT for gfW
    arma::cx_vec FgfW=beta*mcon%arma::fft( mcon%gfWin )%support2;
    // inverse FFT to get fX
    arma::cx_vec Fratio=(FfW/FfU2%support1);
    arma::cx_vec fXF = mcon/beta%arma::ifft( mcon%Fratio);
    // inverse FFT to get gX*fX
    Fratio=(FgfW/FfU2%support1);
    arma::cx_vec gfXF=mcon/beta%arma::ifft( mcon%Fratio);
    // estimate of gX 
    arma::vec ghat = arma::real(gfXF)/arma::real(fXF);
    for (int i=0; i<nx; ++i){
      int ind = (int)(x[i]/beta+m_mid)-1;
      if(ghat[ind]<0.0){
        res(i,jj)=0.0;
      }else{
        res(i,jj) = ghat[ind];
      }
    }
  }
  return List::create(Named("fitxy")=res);
  END_RCPP
}
