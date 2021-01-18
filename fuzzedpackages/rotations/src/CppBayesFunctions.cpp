#include "CppBayesFunctions.h"
#include "basics.h"

Rcpp::NumericVector rcayleyCpp(unsigned int n, double kappa)
{
  Rcpp::RNGScope scope;
  Rcpp::NumericVector bet(n);
  Rcpp::NumericVector alp(n);
  Rcpp::NumericVector theta(n);

  bet = Rcpp::rbeta(n, kappa + 0.5, 1.5);
  alp = Rcpp::rbinom(n, 1, 0.5);

  for (unsigned int i = 0;i < n;++i)
    theta[i] = std::acos(2.0 * bet[i] - 1.0) * (1.0 - 2.0 * alp[i]);

  return theta;
}

double dmbCpp(double r, double kappa)
{
  double den = 2.0 * kappa * std::sqrt(kappa / M_PI);
  den *= std::pow(r, 2.0) * std::exp(-kappa * std::pow(r, 2.0));
  return den;
}

double arsample_mb_unifCpp(double M, double kappa)
{
  Rcpp::RNGScope scope;
  // generate a random observation from target density f assuming g is uniform
  bool found = false;
  Rcpp::NumericVector y(1);
  double evalF = 0.0;
  double x = 0.0;

  while (!found)
  {
    x = Rcpp::as<double>(Rcpp::runif(1, -M_PI, M_PI));
    y = Rcpp::runif(1, 0, M);

    evalF = dmbCpp(x, kappa);

    if (y[0] < evalF)
      found = true;
  }

  return x;
}

Rcpp::NumericVector rar_mb_Cpp(unsigned int n, double kappa, double M)
{
  Rcpp::NumericVector res(n);

  for (unsigned int i = 0;i < n;++i)
    res[i] = arsample_mb_unifCpp(M, kappa);

  return res;
}

Rcpp::NumericVector rmbCpp(unsigned int n, double kappa)
{
  double step = 0.0075;
  double prog = -M_PI;
  double M = 0.0;
  double Mi = 0.0;
  Rcpp::NumericVector res(n);

  while (prog < .5)
  {
    Mi = dmbCpp(prog, kappa);

    if (M < Mi)
      M = Mi;

    prog += step;
  }

  res = rar_mb_Cpp(n, kappa, M);

  return res;
}

double dfisherCpp(double r, double kappa)
{
  if (kappa < 200)
  {
    // Use the matrix Fisher density for kappa<100
    double I02k = R::bessel_i(2.0 * kappa, 0, 1);
    double I12k = R::bessel_i(2.0 * kappa, 1, 1);
    double den = std::exp(2.0 * kappa * std::cos(r));
    den *= (1.0 - std::cos(r));
    den /= (2.0 * M_PI * (I02k - I12k));
    return den;
  }

  // Use the more efficient Maxwell-Boltzman density otherwise
  return dmbCpp(r, kappa);
}

double arsample_unifCpp(double M, double kappa)
{
  Rcpp::RNGScope scope;

  // Generate a random observation from target density f assuming g is uniform
  bool found = false;
  Rcpp::NumericVector y(1);
  double evalF = 0.0;
  double x = 0.0;

  while (!found)
  {
    x = Rcpp::as<double>(Rcpp::runif(1, -M_PI, M_PI));
    y = Rcpp::runif(1, 0, M);

    evalF = dfisherCpp(x, kappa);

    if (y[0] < evalF)
      found = true;
  }

  return x;
}

Rcpp::NumericVector rarCpp(unsigned int n, double kappa, double M)
{
  Rcpp::NumericVector res(n);

  for (unsigned int i = 0;i < n;++i)
    res[i] = arsample_unifCpp(M, kappa);

  return res;
}

Rcpp::NumericVector rfisherCpp(unsigned int n, double kappa)
{
  double step = 0.0075;
  double prog = -M_PI;
  double M = 0.0;
  double Mi = 0.0;

  while (prog < .5)
  {
    Mi = dfisherCpp(prog, kappa);

    if (M < Mi)
      M = Mi;

    prog += step;
  }

  return rarCpp(n, kappa, M);
}

int sign(double x)
{
  if (x < 0)
    return -1;

  return 1;
}

Rcpp::NumericVector rvmisesCPP(unsigned int n, double kappa)
{
  Rcpp::RNGScope scope;
  Rcpp::NumericVector u(3);
  Rcpp::NumericVector theta(n, 10.0);

  u = Rcpp::runif(3, 0, 1);
  double a = 1.0 + std::sqrt(1.0 + 4.0 * std::pow(kappa, 2.0));
  double b = (a - std::sqrt(2.0 * a)) / (2.0 * kappa);
  double r = (1.0 + std::pow(b, 2.0)) / (2.0 * b);
  double z = 0.0;
  double f = 0.0;
  double c = 0.0;

  for (unsigned int i = 0;i < n;++i)
  {
    while (theta[i] > 4)
    {
      // Step 1
      u = Rcpp::runif(3, 0, 1);
      z = std::cos(M_PI * u[0]);
      f = (1.0 + r * z) / (r + z);
      c = kappa * (r - f);

      // Step 2
      u = Rcpp::runif(3, 0, 1);
      if ((c * (2.0 - c) - u[1]) > 0)
        theta[i] = (sign(u[2] - 0.5)) * std::acos(f);
      else
      {
        if ((std::log(c / u[1]) + 1.0 - c) < 0)
          u = Rcpp::runif(3, 0, 1);
        else {
          u = Rcpp::runif(3, 0, 1);
          theta[i] = (sign(u[2] - 0.5)) * std::acos(f);
        }
      }
    }
  }

  return theta;
}

arma::mat centerCpp(const arma::mat &Rs, const arma::mat &S)
{
  // Center the dataset Rs around S, so each row of Rs denoted R is replaces with S'R
  unsigned int n = Rs.n_rows;
  arma::mat Rsi(3, 3);
  arma::mat cRs(n, 9);

  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned int j = 0;j < 9;++j)
     Rsi[j] = Rs(i, j);

    Rsi = S.t() * Rsi;

    for (unsigned int j = 0;j < 9;++j)
      cRs(i, j) = Rsi[j];
  }

  return cRs;
}

double lpvmises(const arma::mat &Rs, const arma::mat &S, double kappa)
{
  // Evaluate the log-posterior for R~von Mises(S,kappa)
  unsigned int n = Rs.n_rows;
  arma::mat cRs = centerCpp(Rs, S);
  arma::mat trcRs(n, 3);

  trcRs.col(0) = cRs.col(0);
  trcRs.col(1) = cRs.col(4);
  trcRs.col(2) = cRs.col(8);

  arma::colvec traces = arma::sum(trcRs, 1);

  double n1 = kappa * (arma::sum(traces) - n) / 2.0;
  double I0k = R::bessel_i(kappa, 0, 1);
  double I1k = R::bessel_i(kappa, 1, 1);
  double n2 = 0.5 * std::log(std::pow(I0k, 2.0) - (I0k * I1k / kappa) - std::pow(I1k, 2.0));
  double d1 = (n + 1.0) * std::log(I0k);
  double d2 = arma::sum(arma::log(3.0 - traces));

  return n1 + n2 - d1 - d2;
}

double lpfisher(const arma::mat &Rs, const arma::mat &S, double kappa)
{
  // Evaluate the log-posterior for R~matrix Fisher(S,kappa)
  unsigned int n = Rs.n_rows;
  arma::mat cRs = centerCpp(Rs, S);
  arma::mat trcRs(n, 3);

  trcRs.col(0) = cRs.col(0);
  trcRs.col(1) = cRs.col(4);
  trcRs.col(2) = cRs.col(8);

  arma::colvec traces = arma::sum(trcRs, 1);

  double n1 = kappa * (arma::sum(traces) - n);
  double I02k = R::bessel_i(2.0 * kappa, 0, 1);
  double I12k = R::bessel_i(2.0 * kappa, 1, 1);
  double n2 = 0.5 * std::log(2.0 * std::pow(I02k, 2.0) / kappa - 2.0 * I02k * I12k / std::pow(kappa, 2.0) + ((1.0 / std::pow(kappa, 2.0)) - (2.0 / kappa)) * std::pow(I12k, 2.0));
  double d1 = (n + 1.0) * std::log(I02k - I12k);

  return n1 + n2 - d1;
}

double lpcayley(const arma::mat &Rs, const arma::mat &S, double kappa)
{
  // Evaluate the log-posterior for R~Cayley(S,kappa)
  unsigned int n = Rs.n_rows;
  arma::mat cRs = centerCpp(Rs, S);
  arma::mat trcRs(n, 3);

  trcRs.col(0) = cRs.col(0);
  trcRs.col(1) = cRs.col(4);
  trcRs.col(2) = cRs.col(8);

  arma::colvec traces = arma::sum(trcRs,  1);

  if (kappa > 169.5)
    kappa = 169.5;

  double p1 = n * std::log(std::sqrt(M_PI) * R::gammafn(kappa + 2.0) / R::gammafn(kappa + 0.5));
  double p2 = 0.5 * std::log(R::trigamma(kappa + 0.5) - R::trigamma(kappa + 2.0));
  double p3 = kappa * arma::sum(arma::log(0.5 + 0.25 * (traces - 1.0)));

  return p1 + p2 + p3;
}

arma::mat genrC(const arma::mat &S, double r)
{
  Rcpp::RNGScope scope;

  Rcpp::NumericVector ta = Rcpp::runif(1, -1, 1);
  double theta = ta[0];
  theta = std::acos(theta);
  double cosr = std::cos(r);
  double sinr = std::sin(r);

  Rcpp::NumericVector ph = Rcpp::runif(1, -M_PI, M_PI);
  double phi = ph[0];

  arma::rowvec u(3);
  arma::colvec ut(3);

  u(0) = std::sin(theta) * std::cos(phi);
  u(1) = std::sin(theta) * std::sin(phi);
  u(2) = std::cos(theta);

  ut = u.t();

  arma::mat Ri;
  arma::mat I(3,3);
  arma::mat SS;
  arma::mat step(3,3);
  I.eye();
  Ri = ut * u;

  step = (I - Ri);
  step *= cosr;
  Ri += step;

  SS = eskewC(u);

  step = SS * sinr;
  Ri += step;
  Ri = S * Ri;

  return Ri;
}

arma::mat S_MCMC_CPP(const arma::mat &Rs,
                     const arma::mat &oldS,
                     double rho,
                     double kappa,
                     int Dist)
{
  Rcpp::RNGScope scope;

  double r, rj1;
  Rcpp::NumericVector W1(1);
  arma::mat Sstar;

  // For now use 'Dist' to identify the form of the likelihood and
  // proposal distribution for S: Cayley is 1, Fisher is 2, von Mises is 3
  // This is ultra inefficient but should work for now

  if (Dist == 1)
  {
    // Generate proposal S~Cayley(oldS,rho) distribution
    r = Rcpp::as<double>(rcayleyCpp(1, rho));
    Sstar = genrC(oldS, r);

    // Compute transition probability: g(Sstar,kappa)/g(oldS,kappa)
    rj1 = lpcayley(Rs, Sstar, kappa);
    rj1 -= lpcayley(Rs, oldS, kappa);
    rj1 = std::exp(rj1);
  }
  else if (Dist == 2)
  {
    // Generate proposal S~Cayley(oldS,rho) distribution
    r = Rcpp::as<double>(rfisherCpp(1, rho));
    Sstar = genrC(oldS,r);

    // Compute transition probability: g(Sstar,kappa)/g(oldS,kappa)
    rj1 = lpfisher(Rs, Sstar, kappa);
    rj1 -= lpfisher(Rs, oldS, kappa);
    rj1 = std::exp(rj1);
  }
  else
  {
    // Generate proposal S~Cayley(oldS,rho) distribution
    r = Rcpp::as<double>(rvmisesCPP(1, rho));
    Sstar = genrC(oldS,r);

    // Compute transition probability: g(Sstar,kappa)/g(oldS,kappa)
    rj1 = lpvmises(Rs, Sstar, kappa);
    rj1 -= lpvmises(Rs, oldS, kappa);
    rj1 = std::exp(rj1);
  }

  if (!std::isfinite(rj1))
    rj1 = 0;

  if (rj1 > 1)
    return Sstar;

  //Generate W~Bern(min(1,rj1)) random variable
  W1 = Rcpp::rbinom(1, 1, rj1);

  if (W1[0] == 1)
    return Sstar;

  return oldS;
}

double kap_MCMC_CPP(const arma::mat &Rs,
                    double oldKappa,
                    double sigma,
                    const arma::mat &S,
                    int Dist)
{
  Rcpp::RNGScope scope;

  double  rj2, kappaStar;
  Rcpp::NumericVector kappaS(1);
  Rcpp::NumericVector W2(1);

  // Generate proposal log(kappa)~N(log(oldKappa),sigma^2)
  kappaS = Rcpp::rnorm(1, std::log(oldKappa), sigma);
  kappaStar = std::exp(kappaS[0]);

  // For now use 'Dist' to identify the form of the likelihood and
  // proposal distribution for S: Cayley is 1, Fisher is 2, von Mises is 3
  // This is ultra inefficient but should work for now

  if (Dist == 1)
  {
    // Compute transition probability:
    rj2 = lpcayley(Rs, S, kappaStar);
    rj2 -= lpcayley(Rs, S, oldKappa);
  }
  else if (Dist == 2)
  {
    // Compute transition probability:
    rj2 = lpfisher(Rs, S, kappaStar);
    rj2 -= lpfisher(Rs, S, oldKappa);
  }
  else
  {
    // Compute transition probability
    rj2 = lpvmises(Rs, S, kappaStar);
    rj2 -= lpvmises(Rs, S, oldKappa);
  }

  rj2 = std::exp(rj2);

  if (!std::isfinite(rj2))
    rj2 = 0;

  if (rj2 > 1)
    return kappaStar;

  //Generate W2~Bern(min(1,rj2))
  W2 = Rcpp::rbinom(1, 1, rj2);

  if (W2[0] == 1)
    return kappaStar;

  return oldKappa;
}

arma::rowvec afun_CPP(const arma::mat &R1, const arma::mat &R2)
{
  unsigned int n = R1.n_rows;
  arma::mat Ri(3, 3);
  arma::rowvec as(n);
  arma::rowvec ds(3);
  as.zeros();

  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned int j = 0;j < 9;++j)
      Ri[j] = R1(i, j);

    Ri = Ri.t() * R2;

    ds(0) = std::acos(Ri(0, 0));
    ds(1) = std::acos(Ri(1, 1));
    ds(2) = std::acos(Ri(2, 2));

    as[i] = arma::max(ds);
  }

  return as;
}

Rcpp::List both_MCMC_CPP(const arma::mat &Rs,
                         arma::mat S0,
                         double kappa0,
                         double rho,
                         double sigma,
                         int burnin,
                         int B,
                         int Dist)
{
  // Rs - the sample
  // S0 - initial central orientation
  // kappa - initial concentration
  // rho - tuning paramter for S proposal distribution, directly related to acceptance rate
  // sigma - tuning paramter for kappa proposal distribution, inversely related to acceptance rate
  // burnin - how may iterations to treat as burnin
  // B - number of draws from posterior to keep
  // Cayley - drawing from Cayley distribution (True) or matrix Fisher (False)
  double Scount = 0.0;
  double Kcount = 0.0;
  arma::mat Sdraws(B, 9);
  Sdraws.zeros();
  Rcpp::NumericVector Kdraws(B);
  arma::mat Snew = S0;
  double Knew = kappa0;
  double Ksame;
  Rcpp::List out;

  for (unsigned int i = 0;i < burnin;++i)
  {
    Snew = S_MCMC_CPP(Rs, Snew, rho, Knew, Dist);
    Knew = kap_MCMC_CPP(Rs, Knew, sigma, Snew, Dist);
  }

  Kdraws[0] = Knew;

  for (unsigned int i = 0;i < 9;++i)
    Sdraws(0, i) = Snew[i];

  for (unsigned int i = 1;i < B;++i)
  {
    S0 = Snew;
    Snew = S_MCMC_CPP(Rs, S0, rho, Kdraws[i - 1], Dist);

    if (arma::accu(arma::abs(S0 - Snew)) < 1e-4)
      Sdraws.row(i) = Sdraws.row(i - 1);
    else
    {
      Scount += 1.0;

      for (unsigned int j = 0;j < 9;++j)
        Sdraws(i, j) = Snew[j];
    }

    Kdraws[i] = kap_MCMC_CPP(Rs, Kdraws[i - 1], sigma, Snew, Dist);
    Ksame = Kdraws[i] - Kdraws[i - 1];

    if (Ksame < 0)
      Ksame *= -1.0;

    if (Ksame > 1e-4)
      Kcount += 1.0;
  }

  out["S"] = Sdraws;
  out["kappa"] = Kdraws;
  out["Saccept"] = Scount / B;
  out["Kaccept"] = Kcount / B;

  return out;
}
