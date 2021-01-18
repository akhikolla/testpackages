/*******************************R-callable function***************************************/
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>  
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

/*============GLOBAL VARIABLES===========*/
/*=====END OF GLOBAL VARIABLES===========*/

/*============FUNCTION DEFINITIONS*/
double mu(double t, double y1, double gamma1, double fH, double f1H, double mu0H, double thetaH, double QH, double *par);
double* func1(double t, double *y, double fH, double f1H, double aH, double bH, double QH, double theta);
void ode45_simpson(double t1, double t2, double y1, double *out, double &s, double nsteps, double fH, double f1H, double aH, double bH, double QH, double thetaH, double mu0H);
double Q(double t, double QH, double theta);
/*=========END OF FUNCTION DEFINITIONS*/

double mu(double t, double y1, double gamma1, double fH, double f1H, double mu0H, double thetaH, double QH) {
  double hfH, hf1H, mu0Ht,  mu;
    
  hfH = fH-y1;
  hf1H = f1H-y1;  
  
  mu0Ht = mu0H*exp(thetaH*t);
  mu = mu0Ht + pow(hfH,2.00)*QH + QH*gamma1;
  
  return mu;
}

//Calculating m (y[1]) & gamma(y[2]):
int func1(double *res, double t, double *y, double fH, double f1H, double aH, double bH, double QH, double theta) {
  double hfH, hf1H, dy1, dy2;
  hfH = fH-y[0];
  hf1H = f1H-y[0];
  dy1 = -1.00*aH*hf1H + 2.00*y[1]*Q(t, QH, theta)*hfH;
  //dy2 = 2.00*aH*y[1] + bH - 2.00*pow(y[1],2.00)*Q(t, QH, theta); //dy2 <- 2*aH*y[2] + bH^2 - 2*y[2]^2*Q(t);
  dy2 = 2.00*aH*y[1] + pow(bH,2) - 2.00*pow(y[1],2.00)*Q(t, QH, theta); //dy2 <- 2*aH*y[2] + bH^2 - 2*y[2]^2*Q(t);
  //double *res = new double[2];
  res[0] = dy1;
  res[1] = dy2;
  
  return 0;
}

double Q(double t, double QH, double theta) {
  double Q;
  Q = QH*exp(theta*t);
  return Q;
}


void ode45_simpson(double t1, double t2, double y1, double *out, double &s, double nsteps, 
          double fH, double f1H, double aH, double bH, double QH, double thetaH, double mu0H) {
  
  double *k1ar, *yfin, *ytmp, *k2ar, *k3ar, *k4ar;
  k1ar = new double[2];
  yfin = new double[2];
  ytmp = new double[2];
  k2ar = new double[2];
  k3ar = new double[2];
  k4ar = new double[2];
  
  double t = t1;
  double h=(t2-t1)/nsteps;
  out[0] = y1;
  out[1] = 0.00;
  
  //Integration:
  s = h/3.00*(-1.00)*mu(t1,y1,0.00, fH, f1H, mu0H, thetaH, QH);
  double ifactor;
      
  for(int j = 0; j < nsteps; j++) {
    //Runge-Kutta method:
    int res = func1(k1ar, t,out, fH, f1H, aH, bH, QH, thetaH);
    yfin[0] = out[0] + h/6.00*k1ar[0];
    yfin[1] = out[1] + h/6.00*k1ar[1];
    ytmp[0] = out[0] + h/2.00*k1ar[0];
    ytmp[1] = out[1] + h/2.00*k1ar[1];
         
    res = func1(k2ar, t,ytmp, fH, f1H, aH, bH, QH, thetaH);
    yfin[0] = yfin[0] + h/3.00*k2ar[0];
    yfin[1] = yfin[1] + h/3.00*k2ar[1];
    ytmp[0] = out[0] + h/2.00*k2ar[0];
    ytmp[1] = out[1] + h/2.00*k2ar[1];
        
    res = func1(k3ar, t,ytmp, fH, f1H, aH, bH, QH, thetaH);
    yfin[0] = yfin[0] + h/3.00*k3ar[0];
    yfin[1] = yfin[1] + h/3.00*k3ar[1];
    ytmp[0] = out[0] + h*k3ar[0];
    ytmp[1] = out[1] + h*k3ar[1];
        
    res = func1(k4ar, t,ytmp, fH, f1H, aH, bH, QH, thetaH);
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
  
  delete k1ar; 
  delete yfin;
  delete ytmp; 
  delete k2ar;
  delete k3ar;
  delete k4ar;  
  
}

RcppExport SEXP complik(SEXP dat, SEXP n, SEXP m, SEXP ah, SEXP f1h, SEXP qh, SEXP bh, SEXP fh, SEXP mu0h, SEXP thetah) {
    
    long N = as<long>(n); //Number of rows
    //long M = as<long>(m); // Number of columns
    double aH = as<double>(ah);
    double f1H = as<double>(f1h);
    double QH = as<double>(qh);
    double bH = as<double>(bh);
    double fH = as<double>(fh);
    double mu0H = as<double>(mu0h); 
    double thetaH  = as<double>(thetah);
    //Actual data set
    Rcpp::NumericMatrix dd = Rcpp::NumericMatrix(dat);   
    
    /*End of data loading*/
    double *out = new double[2];
    double s;
    double  nsteps = 2.00;
    double L; // Likelihood
    L = 0;
    for(int i=0; i<N; i++) {
      //Solving differential equations on intervals:
      double t1 = dd(i,1); 
      double t2 = dd(i,2);
      double y1 = dd(i,3);
      double y2 = dd(i,4);
  
      // Runge-Kutta method:
      ode45_simpson(t1, t2, y1, out, s, nsteps, fH, f1H, aH, bH, QH, thetaH, mu0H);
      
      double m2 = out[0];
      double gamma2 = out[1];
      double pi = 3.141592654;
    
      if(dd(i,0) == 0) { 
        double exp = -0.50*log(2.00*pi*gamma2)-pow((m2-y2),2.00)/2.00/gamma2;
        L = L + s + exp;
        //cout << exp << "\n";
      } else {
        double logprobi = log(1.00 - exp(-1.00*mu(t2, m2, gamma2, fH, f1H, mu0H, thetaH, QH)));
        L = L + s + logprobi;
        //cout << s << " " << logprobi << " " << m2 << " " << gamma2 << "\n";
      }
      //break;
    }
    //std::cout << L << "\n";
    return(Rcpp::wrap(L));
}
