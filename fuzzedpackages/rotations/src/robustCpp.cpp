#include "robustCpp.h"
#include "estimators.h"

arma::rowvec HnCpp(const arma::mat &Qs)
{
  // Compute the Hn tests statistics
  unsigned int n = Qs.n_rows;
  arma::mat T = Qs.t() * Qs;
  arma::mat eigvec;
  arma::mat eigvecJ;
  arma::vec eigval;
  arma::vec eigvalJ;
  arma::eig_sym(eigval, eigvec, T);
  arma::rowvec Hn(n);
  arma::rowvec Qj;
  arma::mat Tj;

  for (unsigned i = 0;i < n;++i)
  {
    Qj = Qs.row(i);
    Tj = T - Qj.t() * Qj;
    arma::eig_sym(eigvalJ, eigvecJ, Tj);
    Hn(i) = (n - 2.0) * (1.0 + eigvalJ(3) - eigval(3)) / (n - 1.0 - eigvalJ(3));
  }

  return Hn;
}

arma::rowvec RdistCArma(const arma::mat &Q1, const arma::rowvec &Q2)
{
  /* Compute the geodesic distance between quaternions Q1 and Q2 */
  /* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single (1x4) quaternion */

	unsigned int n = Q1.n_rows;
	arma::rowvec rs(n);

	for (unsigned int i = 0;i < n;++i)
	{
		double cp = arma::sum(Q1.row(i) * Q2.t());
		rs(i) = std::acos(2.0 * cp * cp - 1.0);
	}

	return rs;
}

arma::rowvec HnCppIntrinsic(const arma::mat &Qs)
{
  // Compute the intrinsic Hn tests statistics

  unsigned int n = Qs.n_rows;

  // Get T matrix of whole sample to make it easier later on
  arma::mat T = Qs.t() * Qs;
  arma::rowvec Qhat = meanQ4C(Qs);
  arma::rowvec dists(n);

  // Sum of squared geometric distances between proj. mean and each obs
  dists = arma::square(RdistCArma(Qs, Qhat));
  double SSE = arma::sum(dists);
  double SSEJ = 0.0;

  arma::rowvec Hn(n);
  arma::rowvec Qhatj;
  arma::mat QsJ(n, 4);

  //Variables for reduced sample mean
  arma::rowvec Qj;
  arma::rowvec distsJ(n - 1);
  arma::mat Tj(4, 4);
  arma::mat eigvecJ(4, 4);

  for (unsigned int i = 0;i < n;++i)
  {
    QsJ.resize(n, 4);

    for (unsigned int j = 0;j < n;++j)
    {
      if (j != i)
        QsJ.row(j) = Qs.row(j);
    }

    QsJ.shed_row(i);

    // Compute projected mean when jth row is cut out
    Qhatj = meanQ4C(QsJ);
    distsJ = arma::square(RdistCArma(QsJ, Qhatj));
    SSEJ = arma::sum(distsJ);

    Hn(i) = (n - 2.0) * (SSE - SSEJ) / SSEJ;
  }

  return Hn;
}

arma::rowvec HnCppBloc(const arma::mat &Qs, const arma::mat &Cs)
{
  // Compute the Hn tests statistics

  unsigned int n = Qs.n_rows;
  unsigned int t = Cs.n_rows;
  unsigned int nc = Cs.n_cols;

  arma::mat T = Qs.t() * Qs;
  arma::mat eigvec;
  arma::mat eigvecJ;
  arma::vec eigval;
  arma::vec eigvalJ;
  arma::eig_sym(eigval, eigvec, T);
  arma::rowvec Hn(nc);
  arma::mat Qj;
  Qj.zeros(t, 4);
  arma::mat::fixed<4,4> Tj;

  for (unsigned int i = 0;i < nc;++i)
  {
    for (unsigned int j = 0;j < t;++j)
    {
      int rowNum = Cs(j, i) - 1;
      Qj.row(j) = Qs.row(rowNum);

      Tj = T - Qj.t() * Qj;
      arma::eig_sym(eigvalJ, eigvecJ, Tj);
      Hn(i) = (n - t - 1.0) * (t + eigvalJ(3) - eigval(3)) / (t * (n - t - eigvalJ(3)));
    }
  }

  return Hn;
}
