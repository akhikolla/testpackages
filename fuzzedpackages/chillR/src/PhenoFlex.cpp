#include <Rcpp.h>
#include <cmath>        
using namespace Rcpp;

const double Pi = 3.14159265358979323846;

inline double P1z(const double T, const double Tu, const double Tb, const double Tc) {
  if(T >= Tb && T <= Tu) {
    return (1./2. * (1 + cos(Pi + Pi * (T - Tb)/(Tu - Tb))) );
  }
  else if(T > Tu && T <= Tc) {
    return ( (1 + cos(Pi/2. + Pi/2. * (T -Tu)/(Tc - Tu))) ); 
  }
  return (0.);
}

inline double P2z(const double T, const double Tu, const double Delta) {
  return( exp(-((T - Tu)/2./Delta)*((T - Tu)/2./Delta)) );
}

inline double PFcn(const double T, const double Tf, const double slope) {
  const double x = slope*Tf*(T-Tf)/T;
  if(x >= 17) return(1);
  else if(x <= -20) return(0);
  const double sr = exp(x);
  return( sr/(1+sr) );
}

//' @title PhenoFlex
//' @description Combined model of the dynamic model for chill accumulation and the GDH model
//'
//' @inheritParams DynModel_driver
//' @param s1 numeric. Slope of transition from chill to heat accumulation
//' @param Tu numeric. GDH optimal temperature 
//' @param Tb numeric. GDH base temperature (lower threshold) 
//' @param Tc numeric. GDH upper temperature (upper threshold)
//' @param yc numeric. Critical value defining end of chill accumulation
//' @param Delta numeric. Width of Gaussian heat accumulation model
//' @param Imodel integer. Heat accumulation model: 0 for GDH and 1 for Gaussian
//' @param zc numeric. Critical value of z determining the end of heat accumulation
//' @param stopatzc boolean. If `TRUE`, the PhenoFlex is applied until the end of the temperature series. Default is to stop once the value zc has been reached.
//' @param basic_output boolean. If `TRUE`, only the bloomindex is returned as a named element of the return list.
//' @useDynLib chillR
//' @author Carsten Urbach <urbach@hiskp.uni-bonn.de>
//' @return
//' A list is returned with named element `bloomindex`, which is the index at which blooming occurs. When `basic_output=FALSE` also `x`, `y`, `z` and `xs` are
//' returned as named element of this list, which are numeric vectors of the same length as the input vector `temp` containing the hourly temperatures.
//' @examples
//' data(KA_weather)
//' hourtemps <- stack_hourly_temps(KA_weather, latitude=50.4)
//' iSeason <- genSeason(hourtemps, years=c(2009))
//' zc <- 190
//' yc <- 40
//' x <- PhenoFlex(temp=hourtemps$hourtemps$Temp[iSeason[[1]]],
//'                times=c(1: length(hourtemps$hourtemps$Temp[iSeason[[1]]])),
//'                zc=zc, stopatzc=TRUE, yc=yc, basic_output=FALSE)
//' DBreakDay <- x$bloomindex
//' ii <- c(1:DBreakDay)
//' plot(x=ii, y=x$z[ii], xlab="Hour Index", ylab="z", col="red", type="l")
//' abline(h=zc, lty=2)
//' plot(x=ii, y=x$y[ii], xlab="Hour Index", ylab="y", col="red", type="l")
//' abline(h=yc, lty=2)
//' @export
// [[Rcpp::export]]
List PhenoFlex(NumericVector temp,
               NumericVector times,
               const double A0=6319.5,
               const double A1=5.939917e13,
               const double E0=3372.8,
               const double E1=9900.3,
               const double slope=1.6,
               const double Tf=4,
               const double s1=0.5,
               const double Tu=25,
               const double Tb=4,
               const double Tc=36,
               const double yc=40,
               const double Delta=4,
               const int Imodel=0,
               const double zc=190,
               bool stopatzc = true,
               bool deg_celsius = true,
               bool basic_output = true) {
  
  const int N = temp.size();
  Rcpp::NumericVector x(N);
  Rcpp::NumericVector y(N);
  Rcpp::NumericVector z(N);
  Rcpp::NumericVector xs(N);
  x[0] = 0.;
  y[0] = 0.;
  z[0] = 0.;
  double _Tf = Tf;
  double _Tu = Tu;
  double _Tc = Tc;
  double _Tb = Tb;
  if(deg_celsius) {
    _Tf += 273.;
    _Tu += 273.;
    _Tc += 273.;
    _Tb += 273.;
  }
  int bloomindex = 0;
  for(int i = 0; i < N-1; i++) {
    double ti = temp[i];
    if(deg_celsius) ti += 273.;
    xs[i] = A0/A1 * exp(-(E0-E1)/ti);
    const double k1 = A1*exp(-E1/ti);
    x[i+1] = xs[i] - (xs[i] - x[i])*exp(-k1*(times[i+1]-times[i]));
    y[i+1] = y[i];
    if(Imodel == 0) {
      z[i+1] = z[i] + P1z(ti, _Tu, _Tb, _Tc) * PFcn(y[i], yc, s1)*(times[i+1]-times[i]);
    }
    else {
      z[i+1] = z[i] + P2z(ti, _Tu, Delta) * PFcn(y[i], yc, s1)*(times[i+1]-times[i]);
    }
    if(x[i+1] >= 1.) {
      double delta = PFcn(ti, _Tf, slope) * x[i+1];
      y[i+1] += delta;
      x[i+1] -= delta;
    }
    if(z[i+1] >= zc) {
      // i+2 for Fortran index convention in R
      bloomindex = i+2;
      if(stopatzc) break;
    }
  }
  if(basic_output) {
    return List::create(Named("bloomindex") = bloomindex);
  }
  return List::create(Named("x") = x, Named("y") = y, Named("z")=z, Named("xs")=xs, Named("bloomindex") = bloomindex);
}
