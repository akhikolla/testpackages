#include "ZhangMethod.h"
#include "estimators.h"

Rcpp::NumericVector RdistC(const Rcpp::NumericMatrix &Q1, const Rcpp::NumericVector &Q2)
{
  /* Compute the geodesic distance between quaternions Q1 and Q2 */
  /* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion */

  unsigned int n = Q1.nrow();
  Rcpp::NumericVector rs(n);

  for (unsigned int i = 0;i < n;++i)
  {
    double cp = Rcpp::sum(Q1(i, Rcpp::_) * Q2);
    rs[i] = std::acos(2.0 * cp * cp - 1.0);
  }

  return rs;
}

arma::rowvec rdistSO3C(const arma::mat &Rs, const arma::mat &R2)
{
  unsigned int n = Rs.n_rows;
  unsigned int m = Rs.n_cols;
  arma::mat R2t = R2.t();

  if (m == 3)
  {
    arma::rowvec theta(1);
    double tri = arma::trace(Rs * R2t);

    if (3.0 - tri < 1.0e-9)
      theta(0) = 0.0;
    else
      theta(0) = std::acos(0.5 * tri - 0.5);

    return theta;
  }

  arma::rowvec theta(n);
  theta.zeros();
  arma::mat33 Rsi;

  for (unsigned int i = 0;i < n;++i)
  {
    for(unsigned int j = 0;j < 9;++j)
      Rsi(j) = Rs(i, j);

    Rsi = Rsi * R2t;
    double tri = arma::trace(Rsi);

    if (3.0 - tri < 1.0e-9)
      theta(i) = 0.0;
    else
      theta(i) = std::acos(0.5 * tri - 0.5);
  }

  return theta;
}

Rcpp::NumericVector EdistC(const Rcpp::NumericMatrix &Q1, const Rcpp::NumericVector &Q2)
{
  /* Compute the Euclidean distance between quaternions Q1 and Q2 */
  /* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion */

  unsigned int n = Q1.nrow();
  Rcpp::NumericVector rs(n);

  for (unsigned int i = 0;i < n;++i)
  {
    double cp = Rcpp::sum(Q1(i, Rcpp::_) * Q2);
    double rsi = 8.0 * (1.0 - cp * cp);
    rs[i] = std::pow(rsi, 0.5);
  }

  return rs;
}

double oneRdistC(const Rcpp::NumericMatrix &Q1, const Rcpp::NumericVector &Q2)
{
  /* Compute the geodesic distance between quaternions Q1 and Q2 */
  /* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion */

  double cp = Rcpp::sum(Q1 * Q2);
  return std::acos(2.0 * cp * cp - 1.0);
}

Rcpp::NumericVector cdfunsC(const Rcpp::NumericMatrix &Qs, const Rcpp::NumericVector &Qhat)
{
  // Compute the projected mean values of c and d to form the pivotal test statistic
  // This estimates c=2E(1-cos(r^2))/3 and d=E(1+2cos(r))/3 from a sample
  unsigned int n = Qs.nrow();

  Rcpp::NumericVector cds(2);
  cds[0] = 0.0;
  cds[1] = 0.0;

  Rcpp::NumericVector rs(n);

  rs = RdistC(Qs, Qhat);

  for (unsigned int i = 0;i < n;++i)
  {
    double crs = std::cos(rs[i]);

    cds[0] += std::pow(crs, 2.0); // c=2E[1-cos(r)^2]/3
    cds[1] += crs;							  // d=E[1+2cos(r)]/3
  }

  cds[0] = 2.0 * (1.0 - cds[0] / n) / 3.0;
  cds[1] = (1.0 + 2.0 * cds[1] / n) / 3.0;

  return cds;
}

Rcpp::NumericVector cdfunsCMedian(const Rcpp::NumericMatrix &Qs,
                                  const Rcpp::NumericVector &Qhat)
{
  // Compute the values c and d to form the pivotal test statistic
  // for the median using quaternions

  unsigned int n = Qs.nrow();

  Rcpp::NumericVector cds(2);
  cds[0] = 0.0;
  cds[1] = 0.0;

  Rcpp::NumericVector rs(n);

  rs = RdistC(Qs, Qhat);

  for (unsigned int i = 0;i < n;++i)
  {
    double crs = std::cos(rs[i]);
    cds[0] += crs; // c=E[1+cos(r)]/6
    double OnemCrs = 1.0 - crs;
    // I think sqrt(1-crs) is close to zero and causing the crash,
    // for now max sure that doesn't happen by taking at least 1e-5
    OnemCrs = std::max(OnemCrs, 1.0e-5);
    cds[1] += (1.0 + 3.0 * crs) * std::pow(OnemCrs, -0.5); // d=E([1+3cos(r)]/12*sqrt[1-cos(r)])
  }

  cds[0] = (1.0 + cds[0] / n) / 6.0;
  cds[1] = cds[1] / n / 12.0;

  return cds;
}

Rcpp::NumericVector zhangQ4(const Rcpp::NumericMatrix &Q, unsigned int m)
{
  unsigned int n = Q.nrow();
  Rcpp::NumericVector cdstar;
  Rcpp::IntegerVector samp(n);
  Rcpp::NumericVector unSamp;

  Rcpp::NumericVector testStat(m);
  arma::mat Qstar(n, 4);
  Rcpp::NumericVector QhatStar;
  Rcpp::NumericMatrix QhatStarMat(1, 4);

  // Convert the sample into armadillo mode
  arma::mat QSamp = Rcpp::as<arma::mat>(Q);

  Rcpp::NumericMatrix QstarRcpp;
  Rcpp::NumericVector Qhat = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(meanQ4C(QSamp)));

  for (unsigned int j = 0;j < m;++j)
  {
    // Bootstrap sample of size n, with replacement
    samp = Rcpp::floor(Rcpp::runif(n, 0, n));
    unSamp = Rcpp::unique(samp);
    unsigned int numUn = unSamp.size();
    unsigned int maxSamp = Rcpp::max(samp);

    while (numUn < 4 || maxSamp > n - 1)
    {
      // If bootstrap samp is less than 4 obs then
      //draw a new sample
      samp = Rcpp::floor(Rcpp::runif(n, 0, n));
      unSamp = Rcpp::unique(samp);
      numUn = unSamp.size();
      maxSamp = Rcpp::max(samp);
    }

    // Copying a matrix row by row produces a bunch of junk messages
    // so I do it with arma instead of standard Rcpp
    for (unsigned int i = 0;i < n;++i)
      Qstar.row(i) = QSamp.row(samp[i]);

    // Both of these functinos return arma variables so
    // They need to be converted to Rcpp type
    QhatStar = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(meanQ4C(Qstar)));
    QstarRcpp = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(Qstar));

    cdstar = cdfunsC(QstarRcpp, QhatStar);
    /* QhatStar needs to be a matrix to be used in RdistC */
    QhatStarMat = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(QhatStar));
    double sqrth = oneRdistC(QhatStarMat, Qhat);

    if (cdstar[0] < 1.0e-7)
      cdstar[0] = 1.0e-7;

    testStat[j] = 2.0 * n * std::pow(cdstar[1], 2.0) * std::pow(sqrth, 2.0) / cdstar[0];
  }

  return testStat;
}

Rcpp::NumericVector cdfunsCSO3(const arma::mat &Rs, const arma::mat &Rhat)
{
  // Compute the projected median values of c and d to form the pivotal test statistic
  // for SO3 data, used by the zhangMedianC function

  unsigned int n = Rs.n_rows;

  Rcpp::NumericVector cds(2);
  cds[0] = 0.0;
  cds[1] = 0.0;

  Rcpp::NumericVector rs(n);
  rs = rdistSO3C(Rs, Rhat);

  for (unsigned int i = 0;i < n;++i)
  {
    double crs = std::cos(rs[i]);

    cds[0] += crs; // c=E[1+cos(r)]/6

    double OnemCrs = 1.0 - crs;
    // I think sqrt(1-crs) is close to zero and causing the crash,
    // for now max sure that doesn't happen by taking at least 1e-5
    OnemCrs = std::max(OnemCrs, 1.0e-5);

    // d=E([1+3cos(r)]/12*sqrt[1-cos(r)])
    cds[1] += (1.0 + 3.0 * crs) * std::pow(OnemCrs, -0.5);
  }

  cds[0] = (1 + cds[0] / n) / 6.0;
  cds[1] = cds[1] / n / 12.0;

  return cds;
}

Rcpp::NumericVector zhangMedianC(const arma::mat &Rs, unsigned int m)
{
  // Compute the bootstrap version of the chang regions for SO3 data
  // because that is what the median function is written for

  // using runif requires this to be set...I think.
  // This has been shown to cause problems in the past so consider
  // using the next line in its place
  Rcpp::RNGScope scope;

  unsigned int n = Rs.n_rows;
  arma::mat Shat = medianSO3C(Rs, 2000, 1.0e-5);
  arma::mat Rstar(n, 9);
  arma::mat Sstar(3, 3);
  Rcpp::NumericVector cdstar(2);
  Rcpp::IntegerVector samp(n);
  Rcpp::NumericVector unSamp;
  Rcpp::NumericVector hsqrtMedian;
  Rcpp::NumericVector hstar(m);

  for (unsigned int j = 0;j < m;++j)
  {
    // Bootstrap sample of size n, with replacement
    samp = Rcpp::floor(Rcpp::runif(n, 0, n));
    unSamp = Rcpp::unique(samp);
    unsigned int numUn = unSamp.size();
    unsigned int maxSamp = Rcpp::max(samp);

    while (numUn < 4 || maxSamp > n - 1)
    {
      // If bootstrap samp is less than 4 obs then
      // draw a new sample
      samp = Rcpp::floor(Rcpp::runif(n, 0, n));
      unSamp = Rcpp::unique(samp);
      numUn = unSamp.size();
      maxSamp = Rcpp::max(samp);
    }

    // Copying a matrix row by row produces a bunch of junk messages
    for (unsigned int i = 0;i < n;++i)
      Rstar.row(i) = Rs.row(samp[i]);

    Sstar = medianSO3C(Rstar, 2000, 1.0e-5);
    cdstar = cdfunsCSO3(Rstar, Sstar);
    hsqrtMedian = rdistSO3C(Shat, Sstar);
    double hsq = hsqrtMedian[0];

    hstar[j] = 2.0 * n * std::pow(cdstar[1], 2.0) * std::pow(hsq, 2.0) / cdstar[0];
  }

  return hstar;
}
