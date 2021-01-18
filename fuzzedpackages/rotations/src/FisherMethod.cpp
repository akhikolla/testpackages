#include "FisherMethod.h"
#include "estimators.h"

double fisherAxisC(const arma::mat &Qs, const arma::rowvec &Qhat)
{
	Rcpp::NumericMatrix Qss = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(Qs));

  unsigned int n = Qs.n_rows;

	arma::mat Qsq = (Qs.t() * Qs) / n;
	arma::mat eigvec;
	arma::mat Mhat(4, 3);
	arma::vec eigval;
	arma::vec eta(3);

  arma::eig_sym(eigval, eigvec, Qsq);

	arma::mat G(3, 3);
	G.zeros();

  for (unsigned int i = 0;i < 3;++i)
  {
    Mhat.col(i) = eigvec.col(i);
    eta(i) = eigval(i);
  }

  for (unsigned int j = 0;j < 3;++j)
  {
		for (unsigned int k = j;k < 3;++k)
		{
			double denom = std::pow(n * (eigval[3] - eta[j]) * (eigval[3] - eta[k]), -1.0);

			for (unsigned int i = 0;i < n;++i)
				G(j, k) = G(j, k) + arma::as_scalar(Qs.row(i) * Mhat.col(j) * Qs.row(i) * Mhat.col(k) * arma::pow(Qs.row(i) * eigvec.col(3), 2.0));

      G(j, k) = G(j, k) * denom;
      G(k, j) = G(j, k);
		}
	}

  arma::mat Ginv = G.i();

  double Tm = arma::as_scalar(n * Qhat * Mhat * Ginv * Mhat.t() * Qhat.t());
  return Tm;
}

double fisherAxisCSymmetric(const arma::mat &Qs, const arma::rowvec &Qhat)
{
	// This is the same as fisherAxisC but the test statistic is much reduced
	// See equation 10 of Fisher et. al. (1996)

	Rcpp::NumericMatrix Qss = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(Qs));

  unsigned int n = Qs.n_rows;

	arma::mat Qsq = (Qs.t() * Qs) / n;
	arma::mat eigvec;
	arma::mat Mhat(4, 3);
	arma::vec eigval;
	arma::vec eta(3);

  arma::eig_sym(eigval, eigvec, Qsq);

  double trGhat = 0.0;
  double Gjj = 0.0;
  double Tnum = 0.0;

  for (unsigned int i = 0;i < 3;++i)
  {
    Mhat.col(i) = eigvec.col(i);
    eta(i) = eigval(i);
  }

  for (unsigned int j = 0;j < 3;++j)
  {
    double denom = std::pow(n * (eigval[3] - eta[j]) * (eigval[3] - eta[j]), -1.0);

    for (unsigned int i = 0;i < n;++i)
      Gjj += arma::as_scalar(Qs.row(i) * Mhat.col(j) * Qs.row(i) * Mhat.col(j) * arma::pow(Qs.row(i) * eigvec.col(3), 2.0));

    Gjj *= denom;
    trGhat += Gjj;
    Gjj = 0.0;
    Tnum += arma::as_scalar(arma::pow(Qhat * Mhat.col(j), 2.0));
  }

  trGhat = std::pow(trGhat, -1.0);

  double Tm = arma::as_scalar(trGhat * (3.0 * n) * Tnum);
  return Tm;
}

arma::vec fisherBootC(const arma::mat &Qs, unsigned int m, bool symm)
{
  unsigned int n = Qs.n_rows;

  arma::rowvec qhat = meanQ4C(Qs);

  arma::vec Tm(m);
  Rcpp::NumericVector unSamp;
  Rcpp::IntegerVector samp(n);
  arma::mat Qstar(n, 4);
  Qstar.zeros();

  for (unsigned int i = 0;i < m;++i)
  {
    samp = Rcpp::floor(Rcpp::runif(n, 0, n));		// Bootstrap sample of size n, with replacement
	  unSamp = Rcpp::unique(samp);
    unsigned int numUn = unSamp.size();

    while (numUn < 4)
    {
      samp = Rcpp::floor(Rcpp::runif(n, 0, n)); // If bootstrap samp is less than 4 obs then
	    unSamp = Rcpp::unique(samp);              // draw a new sample
      numUn = unSamp.size();
    }

    for (unsigned int j = 0;j < n;++j)
      Qstar.row(j) = Qs.row(samp[j]);

    Tm[i] = (symm) ? fisherAxisCSymmetric(Qstar, qhat) : fisherAxisC(Qstar, qhat);
  }

  return Tm;
}
