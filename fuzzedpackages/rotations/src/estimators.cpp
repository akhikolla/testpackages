#include "estimators.h"

int checkQ4(const Rcpp::NumericMatrix &Q)
{
  /* This function will check that the rows in the matrix Q are unit quaternions */
  unsigned int n = Q.nrow();
  unsigned int p = Q.ncol();

  if (n != 4 && p != 4)
  {
    throw Rcpp::exception("The data are not of length 4 each.");
    return 1;
  }

  for (unsigned int i = 0;i < n;++i)
  {
    double len = Rcpp::sum(Q(i, Rcpp::_) * Q(i, Rcpp::_));

    if (len > 1.1 || len < 0.9)
    {
      throw Rcpp::exception("The data are not all unit length so are not quaternions.");
      return 1;
    }
  }

  return 0;
}

int checkSO3(const arma::mat &Rs)
{
  /* This function will check that the rows in the matrix Rs are proper rotations */
  unsigned int n = Rs.n_rows;
  unsigned int p = Rs.n_cols;
  arma::mat Ri(3, 3);
  arma::mat I(3, 3);
  I.eye();

  if (n != 9 && p != 9)
  {
    throw Rcpp::exception("The data are not each of length 9.");
    return 1;
  }

  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned int j = 0;j < 9;++j)
      Ri(j) = Rs(i, j);

    double deti = arma::det(Ri); // Check that each determinant is one

    if (deti > 1.1 || deti < 0.9 )
    {
      throw Rcpp::exception("The data are not all determinant 1, so they are not rotations.");
      return 1;
    }

    double inv = arma::sum(arma::sum(Ri * Ri.t() - I)); // Check that each inverse is the transpose

    if (std::abs(inv) > 0.001)
    {
      throw Rcpp::exception("At least one observation's transpose is not its inverse, so they are not rotations.");
      return 1;
    }
  }

  return 0;
}

arma::mat expskewC(const arma::mat &M)
{
  /* This function takes a 3-by-3 skew symmetric matrix (in so(3)) and
   returns the exponential, a 3-by-3 rotations (in SO(3)) */

  double MMt = arma::sum(arma::sum(M - M.t()));

  if (std::abs(MMt) > 0.01)
    throw Rcpp::exception("The expskewC function is expecting a 3-by-3 skew symmetric matrix.");

  arma::mat expM(3, 3);
  expM.eye();

  double a = std::pow(0.5 * arma::trace(M.t() * M), 0.5);

  if (std::abs(a) < 1.0e-6)
    return expM;

  expM = expM + (std::sin(a) / a) * M + (1.0 - std::cos(a)) * std::pow(a, -2.0) * M * M;

  return expM;
}

arma::mat expskewCMulti(const arma::mat &M)
{
  /* This function takes a sample of 3-by-3 skew symmetric matrices (in so(3)) and
   returns the exponential, a sample of 3-by-3 rotations (in SO(3)) */

  unsigned int n = M.n_rows;
  arma::mat eMi(3, 3);
  arma::mat Mi(3, 3);
  arma::mat expM;
  expM = M; /* take dimensionality of input matrix */
  expM.zeros(); /* make it all zeros */

  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned int j = 0;j < 9;++j)
      Mi(j) = M(i, j);

    eMi = expskewC(Mi);

    for (unsigned int j = 0;j < 9;++j)
      expM(i, j) = eMi(j);
  }

  return expM;
}

arma::mat logSO3C(const arma::mat &R)
{
  arma::mat I(3, 3);
  arma::mat logR(3, 3);
  I.eye();

  double theta = std::acos(0.5 * arma::trace(R) - 0.5);

  if (theta < 0.0001)
  {
    logR.zeros();
    return logR;
  }

  logR = (R - R.t()) * theta / (2.0 * std::sin(theta));

  return logR;
}

arma::mat logSO3CMulti(const arma::mat &R)
{
  // This is a version of logSO3C that allows for multiple rows in Rs

  unsigned int n = R.n_rows;
  arma::mat I(3, 3);
  arma::mat logR(n, 9);
  arma::mat Ri(3, 3);
  arma::mat logRi(3, 3);
  logR.zeros();
  I.eye();
  Ri.zeros();
  logRi.zeros();

  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned j = 0;j < 9;++j)
      Ri(j) = R(i, j);

    double theta = std::acos(0.5 * arma::trace(Ri) - 0.5);

    // If theta<0.0001 leave that row as zeros
    if (theta > 0.00001)
    {
      logRi = (Ri - Ri.t()) * theta / (2.0 * std::sin(theta));

      for (unsigned int j = 0;j < 9;++j)
        logR(i, j) = logRi(j);
    }
  }

  return logR;
}

arma::mat projectSO3C(const arma::mat &M)
{
  /* This function will project the an arbitrary 3-by-3 matrix M in M(3) into SO(3)
   It is expecting a 3-by-3 matrix */

  arma::mat Msq = M.t() * M;
  arma::mat eigvec;
  arma::vec eigval;
  arma::eig_sym(eigval, eigvec, Msq);
  arma::mat dMat(3, 3);
  arma::mat u = arma::fliplr(eigvec);
  dMat.zeros();

  int sign = 1;

  if (arma::det(M) < 0)
    sign = -1;

  dMat(0, 0) = std::pow(eigval[2], -0.5);
  dMat(1, 1) = std::pow(eigval[1], -0.5);
  dMat(2, 2) = sign * std::pow(eigval[0], -0.5);

  return M * u * dMat * u.t();
}

arma::mat meanSO3C(const arma::mat &Rs)
{
  /* Compute the projected mean for a sample of n roations, Rs.
   This function expects Rs to be a n-by-9 matrix where each row
   represents an observations in SO(3) */

  arma::mat Rbarels = arma::mean(Rs);
  arma::mat Rbar(3, 3);

  for (unsigned int i = 0;i < 9;++i)
    Rbar[i] = Rbarels[i];

  return projectSO3C(Rbar);
}

arma::rowvec meanQ4C(const arma::mat &Q)
{
  // Compute the projected mean of the sample Q

  Rcpp::NumericMatrix Qss = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(Q));
  int cq4 = checkQ4(Qss);
  if (cq4)
    throw Rcpp::exception("The data are not in Q4.");

  arma::mat Qsq = Q.t() * Q;
  arma::mat eigvec;
  arma::vec eigval;
  arma::eig_sym(eigval, eigvec, Qsq);
  arma::vec qhat = eigvec.col(3);

  if (qhat[0] < 0)
    qhat = -qhat;

  return qhat.t(); // Want to return it in a row vector so transpose it
}

arma::mat medianSO3C(const arma::mat &Rs, unsigned int maxIterations, double maxEps)
{
  /* Estimate the central direction with the projected median according
   to our algorithm in point estimation paper */

  unsigned int n = Rs.n_rows;
  unsigned int iterations = 0;
  arma::mat S = meanSO3C(Rs);
  arma::mat RsCopy = Rs;
  arma::mat Snew;
  arma::mat33 delta;
  arma::rowvec vnInv(n);
  arma::rowvec deltaV(9);
  arma::rowvec Svec(9);
  double epsilon = 1.0;

  while (epsilon > maxEps && iterations < maxIterations)
  {
    for (unsigned int i = 0;i < 9;++i)
      Svec(i) = S(i);

    double denom = 0;
    for (unsigned int i = 0;i < n;++i)
    {
      vnInv(i) = std::pow(arma::norm(Rs.row(i) - Svec,2), -1.0);
      RsCopy.row(i) = Rs.row(i) * vnInv(i);
      denom += vnInv(i);
    }

    deltaV = arma::sum(RsCopy) / denom;

    for (unsigned int i = 0;i < 9;++i)
      delta(i) = deltaV(i);

    Snew = projectSO3C(delta);

    ++iterations;
    epsilon = arma::norm(Snew - S, 2);
    S = Snew;
  }

  return S;
}

arma::mat HartmedianSO3C(const arma::mat &Rs, unsigned int maxIterations, double maxEps)
{
  /* Estimate the central direction with the projected median according
   to our algorithm in point estimation paper */

  unsigned int n = Rs.n_rows;
  unsigned int iterations = 0;
  arma::mat S = meanSO3C(Rs);
  arma::mat RsCopy = Rs;
  arma::mat Snew;
  arma::mat33 delta;
  arma::mat33 Rsi;
  arma::mat33 vi;
  arma::rowvec vnInv(n);
  double epsilon = 1.0;

  while (epsilon > maxEps && iterations < maxIterations)
  {
    delta.zeros();
    double denom = 0;

    for (unsigned int i = 0;i < n;++i)
    {
      for (unsigned int j = 0;j < 9;++j)
        Rsi(j) = Rs(i, j);

      vi = logSO3C(Rsi * S.t());
      double vin = std::max(arma::norm(vi, 2), 1.0e-5);

      vnInv(i) = std::pow(vin, -1.0);
      delta += vi * vnInv(i);
      denom += vnInv(i);
    }

    delta = delta / denom;
    Snew = expskewC(delta) * S;

    ++iterations;
    epsilon = arma::norm(Snew - S, 2);
    S = Snew;
  }

  return S;
}

arma::mat gmeanSO3C(const arma::mat &Rs, unsigned int maxIterations, double maxEps)
{
  unsigned int n = Rs.n_rows;
  unsigned int iterations = 0;
  arma::mat33 Rsi;
  arma::mat33 r;
  arma::mat S = meanSO3C(Rs);
  double eps = 1.0;
  Rsi.zeros();

  while (eps > maxEps && iterations < maxIterations)
  {
    r.zeros();

    for (unsigned int i = 0;i < n;++i)
    {
      for (unsigned int j = 0;j < 9;++j)
        Rsi(j) = Rs(i, j);

      r = r + logSO3C(S.t() * Rsi);
    }

    r = r / n;
    S = S * expskewC(r);

    eps = arma::norm(r, "fro");
    ++iterations;
  }

  return S;
}
