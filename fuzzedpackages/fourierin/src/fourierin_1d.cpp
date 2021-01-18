
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace arma;

arma::cx_vec fourierin_1d_cpp(const arma::vec & f, double a,
			  double b, double c, double d, double r)
{
  int m = f.n_rows;
  arma::cx_vec out(m), y(2*m), z(2*m), aux(2*m);
  arma::vec J1(m), J2(m), w(m), arg(m);
  double bet, gam, del, cnst;

  bet = (b - a)/m;		// Real numbers
  gam = (d - c)/m;
  del = bet*gam/2;
  J1 = arma::linspace<arma::vec>(0, m-1, m); // mx1 vectors
  J2 = arma::linspace<arma::vec>(m, 2*m-1, m);
  w = c + gam*J1;
  y.zeros();		// (2m) x 1 vector

  // We will first compute the argument and then create the complex
  // vector.
  arg = J1 % (bet*c + del*J1);	// Fill y
  y.rows(0, m - 1) = cx_vec(f % cos(arg), f % sin(arg));

  arg = -del*pow(J1, 2);        // Fill first half of z
  z.rows(0, m - 1) = cx_vec(cos(arg), sin(arg));

  arg = -del*pow(J2 - 2*m, 2);  // Fill 2nd. half of z
  z.rows(m, 2*m - 1) = cx_vec(cos(arg), sin(arg));

  aux = ifft(fft(y) % fft(z));
  out = aux.rows(0, m - 1);
  cnst = bet*pow(2*datum::pi, -(1 - r)/2);
  arg = (a*w + del*pow(J1, 2));
  out = cnst*(out % cx_vec(cos(arg), sin(arg)));

  return out;
}

// [[Rcpp::export]]
arma::cx_vec fourierin_1d_cpp(const arma::vec & f, double a,
			      double b, double c, double d,
			      double r, double s)
{

  int m = f.n_rows;
  arma::cx_vec out(m);

  // fourierin_1d without s argument is meant for s = 1. Thus we have
  // to make it valid for any s.
  out = pow(std::abs(s), 1/2)*fourierin_1d_cpp(f, a, b, s*c, s*d, r);

  return out;
}

// This function compute Fourier integrals evaluated in non-regular
// grids.

// [[Rcpp::export]]
arma::cx_vec fourierin_1d_nonregular_cpp(const arma::vec & f,
					 double a, double b,
					 const arma::vec & w,
					 int resolution,
					 double r, double s)
{
  int m = resolution, k = w.n_rows, i;
  arma::cx_vec out(k);
  arma::vec t(m), arg(m);
  double factor, delta, real, imag;

  delta = (b - a)/m;
  t = arma::linspace<arma::vec>(a + delta/2, b - delta/2, m);
  factor = sqrt(std::abs(s)/pow(2*datum::pi, 1 - r))*delta;

  for(i = 0; i < k; i++)
    {
      arg = s*w(i)*t;
      real = factor*sum(f % cos(arg));
      imag = factor*sum(f % sin(arg));
      out(i) = cx_double(real, imag);
    }

  return out;
}
