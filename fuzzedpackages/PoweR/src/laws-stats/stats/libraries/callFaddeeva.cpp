#include <Rcpp.h>

using namespace Rcpp;

#include "Faddeeva.h"

//' @title Faddeeva family of error functions of the complex variable
//' @description the Faddeeva function
//' @param z complex vector
//' @param relerr double, requested error
//' @return complex vector
//' @describeIn wrap compute w(z) = exp(-z^2) erfc(-iz)
//' @family wrapper
//' @examples 
//' Faddeeva_w(1:10 + 1i)
//' @export
// [[Rcpp::export]]
std::vector< std::complex<double> > Faddeeva_w(const std::vector< std::complex<double> >& z, double relerr=0) { 
  int N = z.size();
  std::vector< std::complex<double> > result(N);
  for(int i=0; i<N; i++) {
    result[i] = Faddeeva::w(z[i], relerr);
  }
  return result;
}

//' the scaled complementary error function
//' @inheritParams Faddeeva_w 
//' @describeIn wrap compute erfcx(z) = exp(z^2) erfc(z)
//' @family wrapper
//' @examples 
//' erfcx(1:10 + 1i)
//' @export
// [[Rcpp::export]]
std::vector< std::complex<double> > erfcx(const std::vector< std::complex<double> >& z, double relerr=0) { 
  int N = z.size();
  std::vector< std::complex<double> > result(N);
  for(int i=0; i<N; i++) {
    result[i] = Faddeeva::erfcx(z[i], relerr);
  }
  return result;
}

//'  the error function of complex arguments
//' @inheritParams Faddeeva_w 
//' @describeIn wrap compute erf(z)
//' @family wrapper
//' @examples 
//' erf(1:10 + 1i)
//' @export
// [[Rcpp::export]]
std::vector< std::complex<double> > erf(const std::vector< std::complex<double> >& z, double relerr=0) { 
  int N = z.size();
  std::vector< std::complex<double> > result(N);
  for(int i=0; i<N; i++) {
    result[i] = Faddeeva::erf(z[i], relerr);
  }
  return result;
}

//' the imaginary error function 
//' @inheritParams Faddeeva_w 
//' @describeIn wrap compute erfi(z) = -i erf(iz)
//' @family wrapper
//' @examples 
//' erfi(1:10 + 1i)
//' @export
// [[Rcpp::export]]
std::vector< std::complex<double> > erfi(const std::vector< std::complex<double> >& z, double relerr=0) { 
  int N = z.size();
  std::vector< std::complex<double> > result(N);
  for(int i=0; i<N; i++) {
    result[i] = Faddeeva::erfi(z[i], relerr);
  }
  return result;
}

//' the complementary error function
//' @inheritParams Faddeeva_w 
//' @describeIn wrap compute erfc(z) = 1 - erf(z)
//' @family wrapper
//' @examples 
//' erfc(1:10 + 1i)
//' @export
// [[Rcpp::export]]
std::vector< std::complex<double> > erfc(const std::vector< std::complex<double> >& z, double relerr=0) { 
  int N = z.size();
  std::vector< std::complex<double> > result(N);
  for(int i=0; i<N; i++) {
    result[i] = Faddeeva::erfc(z[i], relerr);
  }
  return result;
}

//' the Dawson function
//' @inheritParams Faddeeva_w 
//' @describeIn wrap compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
//' @family wrapper
//' @examples 
//' Dawson(1:10 + 1i)
//' @export
// [[Rcpp::export]]
std::vector< std::complex<double> > Dawson(const std::vector< std::complex<double> >& z, double relerr=0) { 
  int N = z.size();
  std::vector< std::complex<double> > result(N);
  for(int i=0; i<N; i++) {
    result[i] = Faddeeva::Dawson(z[i], relerr);
  }
  return result;
}



