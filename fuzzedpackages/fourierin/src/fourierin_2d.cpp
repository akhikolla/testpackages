
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace arma;

arma::cx_mat fourierin_2d_cpp(const arma::mat & f,
			      const arma::vec & a,
			      const arma::vec & b,
			      const arma::vec & c,
			      const arma::vec & d,
			      double r)
{
  int m1 = f.n_rows, m2 = f.n_cols, j1, j2;
  arma::vec m(2), bet(2), gam(2), del(2), a_hat(2),
    w1(m1), w2(m2);
  double arg, aux, cnst;
  arma::cx_vec z1(2*m1), z2(2*m2);
  arma::cx_mat y(2*m1, 2*m2), out(2*m1, 2*m2);

  m[0] = m1;			// Extract dimensions
  m[1] = m2;

  bet = (b - a)/m;		// 2x1 vectors
  gam = (d - c)/m;
  del = bet % gam/2.0;
  a_hat = a + bet/2.0;

  // m1 x 1 vector
  w1 = c[0] + gam[0]*linspace<vec>(0, m1 - 1, m1);

  // m2 x 1 vector
  w2 = c[1] + gam[1]*linspace<vec>(0, m2 - 1, m2);

  // Fill y matrix, m1 x m2
  y.zeros();
  for(j1 = 0; j1 < m1; j1++)
    {
      for(j2 = 0; j2 < m2; j2++)
	{
	arg = j1*(bet[0]*c[0] + del[0]*j1) +
	  j2*(bet[1]*c[1] + del[1]*j2);
	y(j1, j2) = f(j1, j2)*cx_double(cos(arg), sin(arg));
	}
    }

  // Fill z1 (m1 x 1) and z2 (m2 x 1) vectors.
  z1.zeros();			// Initialize vec. with zeros.
  z2.zeros();
  for(j1 = 0; j1 < m1; j1++)	// Start with z1
    {
      j2 = j1 + m1;
      arg = -del(0)*j1*j1;
      z1(j1) = cx_double(cos(arg), sin(arg));
      aux = (j2 - 2.0*m1);
      arg = -1*del(0)*aux*aux;
      z1(j2) = cx_double(cos(arg), sin(arg));
    }
  for(j1 = 0; j1 < m2; j1++)	// Do the same for z2
    {
      j2 = j1 + m2;
      arg = -del(1)*j1*j1;
      z2(j1) = cx_double(cos(arg), sin(arg));
      aux = (j2 - 2.0*m2);
      arg = -del(1)*aux*aux;
      z2(j2) = cx_double(cos(arg), sin(arg));
    }

  // Output m1 x m2 matrix. It only requires to be multiplied by some
  // constants.
  // out = ifft2(fft2(y) % (fft(z1) * (fft(z2).t())));
  out = ifft2(fft2(y) % (fft(z1) * strans(fft(z2))));

  // ... Which we do here.
  cnst = pow(2.0*datum::pi, -(1.0 - r))*
    bet(0)*bet(1);		// Note that n/2 = 1 in this case.
  for(j1 = 0; j1 < m1; j1++)
    {
      for(j2 = 0; j2 < m2; j2++)
  	{
	  arg = (a_hat(0)*w1(j1) + del(0)*j1*j1) + // argument
	    (a_hat(1)*w2(j2) + del(1)*j2*j2);
  	  out(j1, j2) *= cnst*cx_double(cos(arg), sin(arg));
  	}
    }

  return out(span(0, m1 - 1), span(0, m2 - 1));
}

// [[Rcpp::export]]
arma::cx_mat fourierin_2d_cpp(const arma::mat & f,
			      const arma::vec & a,
			      const arma::vec & b,
			      const arma::vec & c,
			      const arma::vec & d,
			      double r, double s)
{

  int m1, m2;
  m1 = f.n_rows;
  m2 = f.n_cols;
  arma::cx_mat out(m1, m2);

  // fourierin_1d without s argument is meant for s = 1. Thus we have
  // to make it valid for any s.
  out = std::abs(s)*fourierin_2d_cpp(f, a, b, s*c, s*d, r);


  return out;
}

/*

  NON-REGULAR GRID INTEGRATION

  To see more details about this function, see documentation in R
  script "fourierin". This function is called by the function
  fourierin_2d.

*/


// This function compute Fourier integrals evaluated in non-regular
// grids.

// [[Rcpp::export]]
arma::cx_mat fourierin_2d_nonregular_cpp(const arma::mat & f,
					 const arma::vec & a,
					 const arma::vec & b,
					 const arma::mat & w,
					 const arma::vec & resolution,
					 double r, double s)
{

  int k = w.n_rows, i, j1, j2;
  arma::cx_vec out(k);
  arma::vec m(resolution), t1(m(0)), t2(m(1)), delta(2);
  double factor, arg;

  delta = (b - a)/m;
  t1 = arma::linspace<arma::vec>(a(0) + delta(0)/2,
				 b(0) - delta(0)/2, m(0));
  t2 = arma::linspace<arma::vec>(a(1) + delta(1)/2,
				 b(1) - delta(1)/2, m(1));
  factor = std::abs(s)/pow(2*datum::pi, 1 - r)*prod(delta);
  out.zeros();

  for(i = 0; i < k; i++)
    {
      for(j1 = 0; j1 < m(0); j1++)
	for(j2 = 0; j2 < m(1); j2++)
	  {
	    arg = s*(t1(j1)*w(i, 0) + t2(j2)*w(i, 1));
	    out(i) += f(j1, j2)*cx_double(cos(arg), sin(arg));
	  }
      out(i) *= factor;
    }

  return out;
}
