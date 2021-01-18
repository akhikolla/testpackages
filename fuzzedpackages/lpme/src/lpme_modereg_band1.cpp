#include "lpme_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// Cross-Validation bandwidth selection without error, based on local constant conditional density
RcppExport SEXP CVdens_LCfit(SEXP X_, SEXP Y_, SEXP pX_, SEXP h1_, SEXP h2_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector pX(pX_);
  NumericVector h1(h1_);
  NumericVector h2(h2_);
  int n = X.size();
  int nh1 = h1.size();
  int nh2 = h2.size();
  
  // temp variable
  NumericMatrix CVh(nh1, nh2);
  NumericVector Ku0XijArray(n*n*nh1);
  arma::cube Ku0Xij(Ku0XijArray.begin(), n, n, nh1, false);
  NumericVector K2YijArray(n*n*nh2);
  arma::cube K2Yij(K2YijArray.begin(), n, n, nh2, false);
  NumericVector intK2YijArray(n*n*nh2);
  arma::cube intK2Yij(intK2YijArray.begin(), n, n, nh2, false);
  
  // Kernel values
  double delta = 0.005; 
  int ngrid1 = round(5.1/delta);
  int ngrid2 = round(7.1/delta);
  NumericVector Ku0x(ngrid1+1);
  NumericVector intK2y(ngrid2+1);
  for(int i=0; i<=ngrid1; ++i){
    Ku0x[i] = exp( -0.5*std::pow(((i+0.0)*delta),2) )/sqrt(2*PI);
  }
  for(int i=0; i<=ngrid2; ++i){
    intK2y[i] = exp( -0.25*std::pow(((i+0.0)*delta),2) )/sqrt(4*PI);
  }
  
  for(int k=0; k<nh1; ++k){
    for(int i=0; i<n; ++i){
      for(int j=i; j<n; ++j){
        int indx1 = round(std::abs(X[i]-X[j])/h1[k]/delta);
        if(indx1>ngrid1) indx1=ngrid1;
        Ku0Xij(i, j, k) = Ku0x[indx1];
        Ku0Xij(j, i, k) = Ku0x[indx1];
      }
    }
  }
  for(int k=0; k<nh2; ++k){
    for(int i=0; i<n; ++i){
      for(int j=i; j<n; ++j){
        int indx1 = round(std::abs(Y[i]-Y[j])/h2[k]/delta);
        int indx2 = indx1;
        if(indx1>ngrid1) indx1=ngrid1;
        if(indx2>ngrid2) indx2=ngrid2;
        K2Yij(i, j, k) = Ku0x[indx1];
        K2Yij(j, i, k) = Ku0x[indx1];
        intK2Yij(i, j, k) = intK2y[indx2];
        intK2Yij(j, i, k) = intK2y[indx2];
      }
    }
  }
  
  for(int k1=0; k1<nh1; ++k1){
    for(int k2=0; k2<nh2; ++k2){
      R_CheckUserInterrupt();
      double cvj=0;
      for(int j=0; j<n; ++j){
        if(pX[j]>1e-10){
          double S0=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              S0 += Ku0Xij(i1, j, k1);
            }
          }
          double num1=0, num2=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              num2 += Ku0Xij(i1, j, k1)*K2Yij(i1,j,k2);
              for(int i2=0; i2<n; ++i2){
                if(i2!=j){
                  num1 += Ku0Xij(i1, j, k1)*Ku0Xij(i2, j, k1)*intK2Yij(i1,i2,k2);
                }
              }
            }
          }
          if(S0<1e-20) S0 = 1e-20;
          cvj += pX[j]*(num1/std::pow(S0,2)-2.0*num2/S0);
        }
      }
      CVh(k1, k2) = cvj/(n+0.0)/h2[k2];
    }
  }
  
  return List::create(Named("CV")=CVh);
  END_RCPP
}

// Cross-Validation bandwidth selection without error, local constant, second order kernel
RcppExport SEXP CVdens_LCfit2(SEXP X_, SEXP Y_, SEXP pX_, SEXP h1_, SEXP h2_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector pX(pX_);
  NumericVector h1(h1_);
  NumericVector h2(h2_);
  int n = X.size();
  int nh1 = h1.size();
  int nh2 = h2.size();
  
  // temp variable
  NumericMatrix CVh(nh1, nh2);
  NumericVector Ku0XijArray(n*n*nh1);
  arma::cube Ku0Xij(Ku0XijArray.begin(), n, n, nh1, false);
  NumericVector K2YijArray(n*n*nh2);
  arma::cube K2Yij(K2YijArray.begin(), n, n, nh2, false);
  NumericVector intK2YijArray(n*n*nh2);
  arma::cube intK2Yij(intK2YijArray.begin(), n, n, nh2, false);
  
  // Kernel values
  double delta = 0.005; 
  int ngrid1 = round(15.0/delta);
  int ngrid2 = round(7.1/delta);
  int ngrid3 = round(5.1/delta);
  NumericVector Ku0x(ngrid1+1);
  NumericVector intK2y(ngrid2+1);
  NumericVector Ku0y(ngrid3+1);
  for(int i=0; i<=ngrid1; ++i){
    Ku0x[i] = K_sec_order((i+0.0)*delta);
  }
  for(int i=0; i<=ngrid2; ++i){
    intK2y[i] = exp( -0.25*std::pow(((i+0.0)*delta),2) )/sqrt(4*PI);
  }
  for(int i=0; i<=ngrid3; ++i){
    Ku0y[i] = exp( -0.5*std::pow(((i+0.0)*delta),2) )/sqrt(2*PI);
  }
  for(int k=0; k<nh1; ++k){
    for(int i=0; i<n; ++i){
      for(int j=i; j<n; ++j){
        int indx1 = round(std::abs(X[i]-X[j])/h1[k]/delta);
        if(indx1>ngrid1) indx1=ngrid1;
        Ku0Xij(i, j, k) = Ku0x[indx1];
        Ku0Xij(j, i, k) = Ku0x[indx1];
      }
    }
  }
  for(int k=0; k<nh2; ++k){
    for(int i=0; i<n; ++i){
      for(int j=i; j<n; ++j){
        int indx1 = round(std::abs(Y[i]-Y[j])/h2[k]/delta);
        int indx2 = indx1;
        if(indx1>ngrid3) indx1=ngrid3;
        if(indx2>ngrid2) indx2=ngrid2;
        K2Yij(i, j, k) = Ku0y[indx1];
        K2Yij(j, i, k) = Ku0y[indx1];
        intK2Yij(i, j, k) = intK2y[indx2];
        intK2Yij(j, i, k) = intK2y[indx2];
      }
    }
  }
  
  for(int k1=0; k1<nh1; ++k1){
    for(int k2=0; k2<nh2; ++k2){
      R_CheckUserInterrupt();
      double cvj=0;
      for(int j=0; j<n; ++j){
        if(pX[j]>1e-10){
          double S0=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              S0 += Ku0Xij(i1, j, k1);
            }
          }
          double num1=0, num2=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              num2 += Ku0Xij(i1, j, k1)*K2Yij(i1,j,k2);
              for(int i2=0; i2<n; ++i2){
                if(i2!=j){
                  num1 += Ku0Xij(i1, j, k1)*Ku0Xij(i2, j, k1)*intK2Yij(i1,i2,k2);
                }
              }
            }
          }
          if(S0<1e-20) S0 = 1e-20;
          cvj += pX[j]*(num1/std::pow(S0,2)-2.0*num2/S0);
        }
      }
      CVh(k1, k2) = cvj/(n+0.0)/h2[k2];
    }
  }
  
  return List::create(Named("CV")=CVh);
  END_RCPP
}

// Cross-Validation bandwidth selection without error, based on local linear conditional density
RcppExport SEXP CVdens_LLfit(SEXP X_, SEXP Y_, SEXP pX_, SEXP h1_, SEXP h2_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector pX(pX_);
  NumericVector h1(h1_);
  NumericVector h2(h2_);
  int n = X.size();
  int nh1 = h1.size();
  int nh2 = h2.size();
  
  // temp variable
  NumericMatrix CVh(nh1, nh2);
  NumericVector Ku0XijArray(n*n*nh1);
  arma::cube Ku0Xij(Ku0XijArray.begin(), n, n, nh1, false);
  NumericVector K2YijArray(n*n*nh2);
  arma::cube K2Yij(K2YijArray.begin(), n, n, nh2, false);
  NumericVector intK2YijArray(n*n*nh2);
  arma::cube intK2Yij(intK2YijArray.begin(), n, n, nh2, false);
  
  // Kernel values
  double delta = 0.005; 
  int ngrid1 = round(5.1/delta);
  int ngrid2 = round(7.1/delta);
  NumericVector Ku0x(ngrid1+1);
  NumericVector intK2y(ngrid2+1);
  for(int i=0; i<=ngrid1; ++i){
    Ku0x[i] = exp( -0.5*std::pow(((i+0.0)*delta),2) )/sqrt(2*PI);
  }
  for(int i=0; i<=ngrid2; ++i){
    intK2y[i] = exp( -0.25*std::pow(((i+0.0)*delta),2) )/sqrt(4*PI);
  }
  
  for(int k=0; k<nh1; ++k){
    for(int i=0; i<n; ++i){
      for(int j=i; j<n; ++j){
        int indx1 = round(std::abs(X[i]-X[j])/h1[k]/delta);
        if(indx1>ngrid1) indx1=ngrid1;
        Ku0Xij(i, j, k) = Ku0x[indx1];
        Ku0Xij(j, i, k) = Ku0x[indx1];
      }
    }
  }
  for(int k=0; k<nh2; ++k){
    for(int i=0; i<n; ++i){
      for(int j=i; j<n; ++j){
        int indx1 = round(std::abs(Y[i]-Y[j])/h2[k]/delta);
        int indx2 = indx1;
        if(indx1>ngrid1) indx1=ngrid1;
        if(indx2>ngrid2) indx2=ngrid2;
        K2Yij(i, j, k) = Ku0x[indx1];
        K2Yij(j, i, k) = Ku0x[indx1];
        intK2Yij(i, j, k) = intK2y[indx2];
        intK2Yij(j, i, k) = intK2y[indx2];
      }
    }
  }
  
  for(int k1=0; k1<nh1; ++k1){
    for(int k2=0; k2<nh2; ++k2){
      R_CheckUserInterrupt();
      double cvj=0;
      for(int j=0; j<n; ++j){
        if(pX[j]>1e-10){
          double S0=0, S1=0, S2=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              double k1xij_tmp = Ku0Xij(i1, j, k1);
              double xij_tmp = ((X[i1]-X[j])/h1[k1]);
              S0 += k1xij_tmp;
              S1 += xij_tmp*k1xij_tmp;
              S2 += std::pow(xij_tmp, 2)*k1xij_tmp;
            }
          }
          double num1=0, num2=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              num2 += (Ku0Xij(i1, j, k1)*S2-((X[i1]-X[j])/h1[k1])*Ku0Xij(i1, j, k1)*S1)*K2Yij(i1,j,k2);
              for(int i2=0; i2<n; ++i2){
                if(i2!=j){
                  num1 += (Ku0Xij(i1, j, k1)*S2-((X[i1]-X[j])/h1[k1])*Ku0Xij(i1, j, k1)*S1)*
                    (Ku0Xij(i2, j, k1)*S2-((X[i2]-X[j])/h1[k1])*Ku0Xij(i2, j, k1)*S1)*intK2Yij(i1,i2,k2);
                }
              }
            }
          }
          double Sden = S0*S2 - std::pow(S1,2);
          if(Sden<1e-20) Sden = 1e-20;
          cvj += pX[j]*(num1/std::pow(Sden,2)-2.0*num2/Sden);
        }
      }
      CVh(k1, k2) = cvj/(n+0.0)/h2[k2];
    }
  }
  
  return List::create(Named("CV")=CVh);
  END_RCPP
}

// Cross-Validation bandwidth selection with laplace error, based on local constant conditional density
RcppExport SEXP CVdens_LCfitLap(SEXP W_, SEXP X_, SEXP Y_, SEXP pX_, SEXP h1_, SEXP h2_,
                                SEXP Ku0x_, SEXP K2y_, SEXP intK2y_, SEXP delta_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector pX(pX_);
  NumericVector h1(h1_);
  NumericVector h2(h2_);
  NumericMatrix Ku0x(Ku0x_);
  NumericVector K2y(K2y_);
  NumericVector intK2y(intK2y_);
  double delta = as<double>(delta_);
  int n = X.size();
  int nh1 = h1.size();
  int nh2 = h2.size();
  int nK1 = Ku0x.nrow()-1;
  int nK2 = K2y.size()-1;
  int nK3 = intK2y.size()-1;
  
  // temp variable
  NumericMatrix CVh(nh1, nh2);
  NumericVector Ku0XijArray(n*n*nh1);
  arma::cube Ku0Xij(Ku0XijArray.begin(), n, n, nh1, false);
  NumericVector K2YijArray(n*n*nh2);
  arma::cube K2Yij(K2YijArray.begin(), n, n, nh2, false);
  NumericVector intK2YijArray(n*n*nh2);
  arma::cube intK2Yij(intK2YijArray.begin(), n, n, nh2, false);
  
  for(int k=0; k<nh1; ++k){
    for(int i=0; i<n; ++i){
      for(int j=0; j<n; ++j){
        int indx = round(std::abs(W[i]-X[j])/h1[k]/delta);
        if(indx>nK1) indx=nK1;
        Ku0Xij(i, j, k) = Ku0x(indx,k);
      }
    }
  }
  for(int k=0; k<nh2; ++k){
    for(int i=0; i<n; ++i){
      for(int j=i; j<n; ++j){
        int indx2 = round(std::abs(Y[i]-Y[j])/h2[k]/delta);
        int indx3 = indx2;
        if(indx2>nK2) indx2=nK2;
        if(indx3>nK3) indx3=nK3;
        K2Yij(i, j, k) = K2y[indx2];
        K2Yij(j, i, k) = K2y[indx2];
        intK2Yij(i, j, k) = intK2y[indx3];
        intK2Yij(j, i, k) = intK2y[indx3];
      }
    }
  }
  
  for(int k1=0; k1<nh1; ++k1){
    for(int k2=0; k2<nh2; ++k2){
      R_CheckUserInterrupt();
      double cvj=0;
      for(int j=0; j<n; ++j){
        if(pX[j]>1e-10){
          double S0=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              S0 += Ku0Xij(i1, j, k1);
            }
          }
          double num1=0, num2=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              num2 += Ku0Xij(i1, j, k1)*K2Yij(i1,j,k2);
              for(int i2=0; i2<n; ++i2){
                if(i2!=j){
                  num1 += Ku0Xij(i1, j, k1)*Ku0Xij(i2, j, k1)*intK2Yij(i1,i2,k2);
                }
              }
            }
          }
          if(S0<1e-20) S0 = 1e-20;
          cvj += pX[j]*(num1/std::pow(S0,2)-2.0*num2/S0);
        }
      }
      CVh(k1, k2) = cvj/(n+0.0)/h2[k2];
    }
  }
  
  return List::create(Named("CV")=CVh);
  END_RCPP
}

// Cross-Validation bandwidth selection with laplace error, based on local linear conditional density
RcppExport SEXP CVdens_LLfitLap(SEXP W_, SEXP X_, SEXP Y_, SEXP pX_, SEXP h1_, SEXP h2_,
                                SEXP Ku0x_, SEXP Ku1x_, SEXP Ku2x_, SEXP K2y_, SEXP intK2y_, SEXP delta_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  NumericVector W(W_);
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector pX(pX_);
  NumericVector h1(h1_);
  NumericVector h2(h2_);
  NumericMatrix Ku0x(Ku0x_);
  NumericMatrix Ku1x(Ku1x_);
  NumericMatrix Ku2x(Ku2x_);
  NumericVector K2y(K2y_);
  NumericVector intK2y(intK2y_);
  double delta = as<double>(delta_);
  int n = X.size();
  int nh1 = h1.size();
  int nh2 = h2.size();
  int nK1 = Ku0x.nrow()-1;
  int nK2 = K2y.size()-1;
  int nK3 = intK2y.size()-1;
  
  // temp variable
  NumericMatrix CVh(nh1, nh2);
  NumericVector Ku0XijArray(n*n*nh1);
  arma::cube Ku0Xij(Ku0XijArray.begin(), n, n, nh1, false);
  NumericVector Ku1XijArray(n*n*nh1);
  arma::cube Ku1Xij(Ku1XijArray.begin(), n, n, nh1, false);
  NumericVector Ku2XijArray(n*n*nh1);
  arma::cube Ku2Xij(Ku2XijArray.begin(), n, n, nh1, false);
  NumericVector K2YijArray(n*n*nh2);
  arma::cube K2Yij(K2YijArray.begin(), n, n, nh2, false);
  NumericVector intK2YijArray(n*n*nh2);
  arma::cube intK2Yij(intK2YijArray.begin(), n, n, nh2, false);
  
  for(int k=0; k<nh1; ++k){
    for(int i=0; i<n; ++i){
      for(int j=0; j<n; ++j){
        int indx = round(std::abs(W[i]-X[j])/h1[k]/delta);
        if(indx>nK1) indx=nK1;
        Ku0Xij(i, j, k) = Ku0x(indx,k);
        if((W[i]-X[j])<0){
          Ku1Xij(i, j, k) = 0.0-Ku1x(indx,k);
        }else{
          Ku1Xij(i, j, k) = Ku1x(indx,k);
        }
        Ku2Xij(i, j, k) = Ku2x(indx,k);
      }
    }
  }
  for(int k=0; k<nh2; ++k){
    for(int i=0; i<n; ++i){
      for(int j=i; j<n; ++j){
        int indx2 = round(std::abs(Y[i]-Y[j])/h2[k]/delta);
        int indx3 = indx2;
        if(indx2>nK2) indx2=nK2;
        if(indx3>nK3) indx3=nK3;
        K2Yij(i, j, k) = K2y[indx2];
        K2Yij(j, i, k) = K2y[indx2];
        intK2Yij(i, j, k) = intK2y[indx3];
        intK2Yij(j, i, k) = intK2y[indx3];
      }
    }
  }
  
  for(int k1=0; k1<nh1; ++k1){
    for(int k2=0; k2<nh2; ++k2){
      R_CheckUserInterrupt();
      double cvj=0;
      for(int j=0; j<n; ++j){
        if(pX[j]>1e-10){
          double S0=0, S1=0, S2=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              S0 += Ku0Xij(i1, j, k1);
              S1 += Ku1Xij(i1, j, k1);
              S2 += Ku2Xij(i1, j, k1);
            }
          }
          double num1=0, num2=0;
          for( int i1=0; i1<n; ++i1){
            if(i1!=j){
              num2 += (Ku0Xij(i1, j, k1)*S2-Ku1Xij(i1, j, k1)*S1)*K2Yij(i1,j,k2);
              for(int i2=0; i2<n; ++i2){
                if(i2!=j){
                  num1 += (Ku0Xij(i1, j, k1)*S2-Ku1Xij(i1, j, k1)*S1)*
                    (Ku0Xij(i2, j, k1)*S2-Ku1Xij(i2, j, k1)*S1)*intK2Yij(i1,i2,k2);
                }
              }
            }
          }
          double Sden = S0*S2 - std::pow(S1,2);
          if(Sden<1e-20) Sden = 1e-20;
          cvj += pX[j]*(num1/std::pow(Sden,2)-2.0*num2/Sden);
        }
      }
      CVh(k1, k2) = cvj/(n+0.0)/h2[k2];
    }
  }
  
  return List::create(Named("CV")=CVh);
  END_RCPP
}
