#include "basics.h"

arma::mat eskewC(const arma::rowvec &U)
{
  double ulen = arma::norm(U, 2);
  double u = U(0) / ulen;
  double v = U(1) / ulen;
  double w = U(2) / ulen;

  arma::mat res(3,3);
  res.zeros();
  res(0, 1) = -w;
  res(1, 0) =  w;
  res(0, 2) =  v;
  res(2, 0) = -v;
  res(1, 2) = -u;
  res(2, 1) =  u;

  return res;
}

arma::mat SO3defaultC(const arma::mat &U, const arma::vec &theta)
{
  // U is an n-by-3 matrix, each row is a misorentation axis
  // theta is a vector of length n, each item is a misorientation angle
  unsigned int n = U.n_rows;
  arma::mat Ri(3, 3);
  arma::mat Rs(n, 9);
  arma::mat I(3, 3);
  arma::mat SS(3, 3);
  arma::rowvec Rir;
  Ri.zeros();
  Rs.zeros();
  I.eye();

  for (unsigned int i = 0;i < n;++i)
  {
    Ri = U.row(i).t() * U.row(i);
    SS = eskewC(U.row(i));
    Ri = Ri + (I - Ri) * std::cos(theta[i]) + SS * std::sin(theta[i]);
    Rs.row(i) = Rcpp::as<arma::rowvec>(Rcpp::wrap(Ri));
  }

  return Rs;
}

arma::mat Q4defaultC(const arma::mat &U, const arma::vec &theta)
{
  unsigned int n1 = U.n_rows;
  unsigned int n = theta.size();
  arma::mat q(n, 4);
  q.zeros();

  if (n1 != n)
    return q;

  arma::vec stheta = arma::sin(theta / 2.0);

  q.col(0) = arma::cos(theta / 2.0);
  q.col(1) = U.col(0) % stheta;
  q.col(2) = U.col(1) % stheta;
  q.col(3) = U.col(2) % stheta;

  return q;
}

arma::mat pMatC(const arma::mat &p)
{
  arma::mat Pmat(4, 4);
  arma::mat revI(4, 4);
  Pmat.zeros();
  revI.zeros();

  Pmat.col(0) = arma::conv_to<arma::vec>::from(p);

  revI(0, 1) = -1;
  revI(1, 0) =  1;
  revI(2, 3) =  1;
  revI(3, 2) = -1;
  Pmat.col(1) = revI * p;

  revI.zeros();
  revI(0, 2) = -1;
  revI(1, 3) = -1;
  revI(2, 0) =  1;
  revI(3, 1) =  1;
  Pmat.col(2) = revI * p;

  revI.zeros();
  revI(0, 3) = -1;
  revI(1, 2) =  1;
  revI(2, 1) = -1;
  revI(3, 0) =  1;
  Pmat.col(3) = revI * p;

  return Pmat;
}

arma::mat genrC(const arma::vec &r, const arma::mat &S, int SO3, const arma::mat &u)
{
  // r is a vector of angles
  // S is the central direction
  // SO3 is an integer, 1 means SO3, anything else gives
  Rcpp::RNGScope scope;

  unsigned int n = r.size();
  unsigned int n1 = u.n_rows;
  unsigned int n2 = u.n_cols;

  if (n1 != n || n2 != 3)
  {
    arma::mat q(n, 4);
    q.zeros();
    return q;
  }

  if (SO3 == 1)
  {
    arma::mat Rs(n, 9);
    arma::mat33 Rsi;
    Rs = SO3defaultC(u, r);

    for (unsigned int i = 0;i < n;++i)
    {
      for (unsigned int j = 0;j < 9;++j)
        Rsi(j) = Rs(i, j);

      Rsi = S * Rsi;
      Rs.row(i) = Rcpp::as<arma::rowvec>(Rcpp::wrap(Rsi));
    }

    return Rs;
  }

  arma::mat q;
  q.zeros(n, 4);
  q = Q4defaultC(u, r);
  return q;
}
