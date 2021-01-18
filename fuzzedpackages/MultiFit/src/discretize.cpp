#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::List discretizeCpp( arma::mat a,
                          arma::mat b,
                          arma::rowvec w,
                          arma::uvec mask,
                          arma::umat ij,
                          int Dx,
                          int Dy
                       ) {
  // # Examine only data points from parent cuboid:
  arma::uvec msk = find(mask>0);
  a = a.rows(msk);
  b = b.rows(msk);

  int m = a.n_rows;
  int o = accu(mask);

  rowvec lx = w.subvec( Dx+Dy, Dx+Dy+Dx-1 );
  rowvec kx = w.subvec( 0, Dx-1 );
  rowvec p2kx(Dx);
  for (int i=0; i<Dx; i++) p2kx(i) = pow(2, -kx(i));
  rowvec xl = (lx-1.0) % p2kx;
  rowvec xh = lx % p2kx;

  rowvec ly = w.subvec( Dx+Dy+Dx, Dx+Dy+Dx+Dy-1 );
  rowvec ky = w.subvec( Dx, Dx+Dy-1 );
  rowvec p2ky(Dy);
  for (int j=0; j<Dy; j++) p2ky(j) = pow(2, -ky(j));
  rowvec yl = (ly-1.0) % p2ky;
  rowvec yh = ly % p2ky;

  // # For each pair of margins, discretize to a 2x2 contingecy table:
  rowvec xm = (xl+xh)/2.0;
  mat x_mid = repmat(xm, m, 1);
  rowvec ym = (yl+yh)/2.0;
  mat y_mid = repmat(ym, m, 1);
  umat x0 = a < x_mid ;
  umat y0 = b < y_mid ;
  umat tables(ij.n_rows, 8);
  tables.cols(0, 1) = ij;
  for (unsigned int c=0; c<ij.n_rows; c++) {
    int i = ij(c, 0) - 1;
    int j = ij(c, 1) - 1;
    tables(c,2) = accu(x0.col(i) % y0.col(j));
    tables(c,3) = accu(x0.col(i) % (1-y0.col(j)));
    tables(c,4) = accu((1-x0.col(i)) % y0.col(j));
    tables(c,5) = o - tables(c,2) - tables(c,3) - tables(c,4);
  }
  tables.col(6).fill(w(w.n_cols - 1));
  tables.col(7).fill(0);

  return Rcpp::List::create(
    Rcpp::Named( "tables" ) = tables,
    Rcpp::Named( "mask" ) = mask,
    Rcpp::Named( "x0.sub.mask" ) = x0,
    Rcpp::Named( "y0.sub.mask" ) = y0 );
}
