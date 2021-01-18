#include <iostream>
#include "RcppArmadillo.h"
#include "../inst/include/Countr_types.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Paper refers to the The Paper by Baker & Kharrat: Discrete
// distributions from renewal processes: fast computation of probabilities
// This file corresponds to Section 2

// get the convolution of 2 pdfs of different order
//
// convolves the pdf of some number of events with another pdf to get the new
// order pdf.
//
// This routine convolves \code{p}, the pdf of some number
// of events occurring, with \code{df} (ditto) to get the resulting
// new order pdf. See section 2 of the paper for more details.
// @param p double vector (passed by reference in the c++ code). First pdf
// @param df double vector (passed by cte reference in the c++ code)
// @param nprobs unsigned integer maximum count required. All probabilities
// from 0 up to \code{nprobs} will be returned.
// @param nsteps integer number of steps used.
// @param probs double vector vector of probabilities to be computed.
// This vector will be passed by reference and will be updated.
// @return only \code{p} and \code{probs} are updated and nothing is returned.
arma::vec convolve(unsigned nprobs, const arma::vec& df, arma::vec& p,
		   const unsigned& nsteps) {
  unsigned klow = 1;
  unsigned n, np, k;
  double ptemp;
  arma::vec probs(nprobs + 1, fill::zeros);

  for(n = 0; n < nprobs; n++) {
    np = n + 1;
    if (np == nprobs)
      klow = nsteps;

    for(k = nsteps; k >= klow; k --) {
      ptemp = 0.0;
      for (unsigned j = 1; j<=k; j++)
	ptemp = ptemp + p(k - j + 1) * df(j);

      p(k) = ptemp;
    }
    probs(np) = p(nsteps);
    if (np != nprobs) {
      for(k = nsteps; k >= 1; k --)
	p(k) = 0.5 * (p(k) + p(k - 1));
    }
  }

  return(probs);
}

arma::vec doOneConvolution(unsigned xmax, arma::vec& p,
			   arma::vec& df,
			   const arma::vec& fwork, const unsigned& nsteps,
			   const unsigned fact) {
  double sth;
  unsigned i, ik;
  double stl = 1.0;
  unsigned fact2 = fact * 0.5;
  for(i = 1 ;i <= nsteps; i++) {
      ik = fact * i;
      sth = fwork(ik);
      df(i) = stl - sth;
      stl = sth;
      p(i) = fwork(ik - fact2);
  }

  return(convolve(xmax, df, p, nsteps));
}

// Compute all probabilities up to \code{xmax} by convolution
//
// Compute all probabilities up to \code{xmax} by convolution (and eventually
// improved by Richardson extrapolation).
//
// The routine does convolutions to produce probabilities \code{probs(0)},
// ... \code{probs(xmax)} using \code{nsteps} steps, and refines result by
// Richardson extrapolation if \code{extrap} is \code{TRUE}.
// @param xmax unsigned integer maximun probability required.
// @param distPars Rcpp::List list of parameters for the desired distribution
// @param extrapolPars arma::vec of length 2. The extrapolation values.
// @param nsteps unsiged integer number of steps used to compute the integral.
// @param time double time at wich to compute the probabilities. Set to 1 by
// default.
// @param extrap logical if \code{TRUE}, Richardson extrapolation will be
// applied to improve accuracy.
// @return vector of probabilities \code{probs(0)}, ... \code{probs(xmax)}.
arma::vec getAllProbs(unsigned xmax, const Rcpp::List distPars,
		      arma::vec extrapolPars, const std::string dist,
		      const unsigned& nsteps = 100,
		      double time = 1.0, bool extrap = true) {

  arma::vec probs(xmax + 1, fill::zeros);

  double stl = 1.0;
  double sth = 1.0;
  double en = (double) nsteps;
  // unused: double h = time / en;
  double xi, th, tee;
  unsigned i;
  vec df, p;

  if (extrap) { // use Richardson extrapolation to reduce the error
    // define the steps needed
    unsigned nsteps1 = nsteps / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    p = zeros<vec> (needed + 1);
    df = zeros<vec> (needed + 1);
    vec fwork(needed + 1, fill::zeros);

    for(i = 1; i <= needed; i++) {
      xi = double (i);
      tee = time * xi / en;
      fwork(i) = surv(tee, distPars, dist);
    }

    // =========================== concolutions ================================
    // ---------- conv1
    arma::vec probs1 = doOneConvolution(xmax, p, df, fwork, nsteps1, 8);
    // ---------- conv2
    arma::vec probs2 = doOneConvolution(xmax, p, df, fwork, nsteps2, 4);
    // ---------- conv3
    probs = doOneConvolution(xmax, p, df, fwork, nsteps3, 2);
    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
    probs(0) = fwork(needed);
  } else {
    p = zeros<vec> (nsteps + 1);
    df = zeros<vec> (nsteps + 1);

    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      th = time * xi / en;
      sth = surv(th, distPars, dist);;
      df(i) = (stl - sth);
      p(i) = surv(th - 0.5 * time / en, distPars, dist);
      stl = sth;
    }
    probs = convolve(xmax, df, p, nsteps);
    probs(0) = sth;
  }

  return(probs);
}

arma::vec getAllProbs(unsigned xmax, const Rcpp::List distPars,
		      arma::vec extrapolPars, Rcpp::Function survR,
		      const unsigned& nsteps = 100,
		      double time = 1.0, bool extrap = true) {

  arma::vec probs(xmax + 1, fill::zeros);

  double stl = 1.0;
  double sth = 1.0;
  double en = (double) nsteps;
  // double h = time / en;
  double xi, th, tee;
  unsigned i;
  vec df, p;
  Rcpp::NumericVector rTemp;

  if (extrap) { // use Richardson extrapolation to reduce the error
    // unsigned i8, i4, i2;
    // define the steps needed
    unsigned nsteps1 = nsteps / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    p = zeros<vec> (needed + 1);
    df = zeros<vec> (needed + 1);
    vec fwork(needed + 1, fill::zeros);

    for(i = 1; i <= needed; i++) {
      xi = double (i);
      tee = time * xi / en;
      rTemp = survR(tee, distPars);
      fwork(i) = rTemp[0];
    }

    // =========================== concolutions ================================
    // ---------- conv1
    arma::vec probs1 = doOneConvolution(xmax, p, df, fwork, nsteps1, 8);
    // ---------- conv2
    arma::vec probs2 = doOneConvolution(xmax, p, df, fwork, nsteps2, 4);
    // ---------- conv3
    probs = doOneConvolution(xmax, p, df, fwork, nsteps3, 2);
    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
    probs(0) = fwork(needed);
  } else {
    p = zeros<vec> (nsteps + 1);
    df = zeros<vec> (nsteps + 1);

    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      th = time * xi / en;
      rTemp = survR(th, distPars);
      sth = rTemp[0];
      df(i) = (stl - sth);
      rTemp = survR(th - 0.5 * time / en, distPars);
      p(i) = rTemp[0];
      stl = sth;
    }
    probs = convolve(xmax, df, p, nsteps);
    probs(0) = sth;
  }

  return(probs);
}

//' Compute count probabilities using simple convolution
//'
//' Compute count probabilities using simple convolution (section 2) for the
//' built-in distributions
//'
//' The routine does convolutions to produce probabilities \code{probs(0)},
//' ... \code{probs(xmax)} using \code{nsteps} steps, and refines result by
//' Richardson extrapolation if \code{extrap} is \code{TRUE} using the
//' algorithm of section 2.
//'
//' @param x integer (vector), the desired count values.
//' @inheritParams surv
//' @param nsteps unsiged integer number of steps used to compute the integral.
//' @param time double time at wich to compute the probabilities. Set to 1 by
//' default.
//' @param extrap logical if \code{TRUE}, Richardson extrapolation will be
//' applied to improve accuracy.
//' @param logFlag logical if \code{TRUE} the log-probability will be returned.
//' @return vector of probabilities P(x(i)) for i = 1, ..., n where n is
//' \code{length} of \code{x}.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_allProbs_bi(arma::Col<unsigned> x, const Rcpp::List distPars,
			     const std::string dist,
			     const unsigned& nsteps = 100,
			     double time = 1.0, bool extrap = true,
			     bool logFlag = false) {

  double xmax = max(x);
  arma::vec extrapolPars = getextrapolPars(distPars, dist);
  arma::vec all = getAllProbs(xmax, distPars, extrapolPars, dist, nsteps,
			      time, extrap);
  arma::vec vals(x.n_elem, fill::zeros);
  arma::Col<unsigned> x_unique = unique(x);
  double xk;

   for (unsigned k = 0; k < x_unique.n_elem; k ++) {
     xk = x_unique(k);
     uvec ind = find(x == xk);
     arma::vec valsk(ind.n_elem, fill::ones);
     vals.elem(ind) = valsk * all(xk);
   }

    if (logFlag)
    return(log(vals));
  else
    return(vals);
}

//' Compute count probabilities using simple convolution
//'
//' Compute count probabilities using simple convolution (section 2) for user
//' passed survival functions
//'
//' @param extrapolPars ma::vec of length 2. The extrapolation values.
//' @param survR Rcpp::Function user passed survival function; should have the
//' signature \code{function(t, distPars)} where \code{t} is a real number (>0)
//' where the survival function is evaluated and \code{distPars} is a list of
//' distribution parameters. It should return a double value.
//' @inheritParams dCount_allProbs_bi
//' @rdname dCount_allProbs_bi
//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_allProbs_user(arma::Col<unsigned> x, const Rcpp::List distPars,
			      arma::vec extrapolPars, Rcpp::Function survR,
			      const unsigned& nsteps = 100,
			      double time = 1.0, bool extrap = true,
			      bool logFlag = false) {

  double xmax = max(x);
  arma::vec all = getAllProbs(xmax, distPars, extrapolPars, survR, nsteps,
			      time, extrap);
  arma::vec vals(x.n_elem, fill::zeros);
  arma::Col<unsigned> x_unique = unique(x);
  double xk;

   for (unsigned k = 0; k < x_unique.n_elem; k ++) {
     xk = x_unique(k);
     uvec ind = find(x == xk);
     arma::vec valsk(ind.n_elem, fill::ones);
     vals.elem(ind) = valsk * all(xk);
   }

    if (logFlag)
    return(log(vals));
  else
    return(vals);
}

//' @keywords internal
// [[Rcpp::export]]
double dCount_allProbs_scalar_bi(unsigned x, const Rcpp::List distPars,
				 const std::string dist,
				 const unsigned& nsteps = 100,
				 double time = 1.0, bool extrap = true,
				 bool logFlag = false) {

  arma::vec extrapolPars = getextrapolPars(distPars, dist);
  arma::vec all = getAllProbs(x, distPars, extrapolPars, dist, nsteps,
			      time, extrap);

  double out = all(all.n_elem - 1);
  if (logFlag)
    out = log(out);

  return(out);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_allProbs_vec_bi(arma::Col<unsigned> x, const Rcpp::List distPars,
				 const std::string dist,
				 const unsigned& nsteps = 100,
				 double time = 1.0, bool extrap = true,
				 bool logFlag = false) {
  // 2018-04-12 was: unsigned lnt = x.n_elem;
  int lnt = x.n_elem;
  arma::vec pbs(lnt, fill::zeros);
  Rcpp::List distParsi;

  if (lnt != distPars.size())
    stop("x and distPars should have same length !");

  for (int i = 0; i < lnt; i++) {
    distParsi = distPars[i];
    pbs[i] = dCount_allProbs_scalar_bi(x[i], distParsi, dist, nsteps, time,
				       extrap, logFlag);
  }

  return(pbs);
}

//' @keywords internal
// [[Rcpp::export]]
double dCount_allProbs_scalar_user(unsigned x, const Rcpp::List distPars,
				  arma::vec extrapolPars, Rcpp::Function survR,
				  const unsigned& nsteps = 100,
				  double time = 1.0, bool extrap = true,
				  bool logFlag = false) {

  arma::vec all = getAllProbs(x, distPars, extrapolPars, survR, nsteps,
			      time, extrap);

  double out = all(all.n_elem - 1);
  if (logFlag)
    out = log(out);

  return(out);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_allProbs_vec_user(arma::Col<unsigned> x,
				   const Rcpp::List distPars,
				   const Rcpp::List extrapolPars,
				   Rcpp::Function survR,
				   const unsigned& nsteps = 100,
				   double time = 1.0, bool extrap = true,
				   bool logFlag = false) {
  // 2018-04-12 was: unsigned lnt = x.n_elem;
  int lnt = x.n_elem;
  arma::vec pbs(lnt, fill::zeros);
  Rcpp::List distParsi;

  if (lnt != distPars.size())
    stop("x and distPars should have same length !");

  if (lnt != extrapolPars.size())
    stop("x and distPars should have same length !");

  for (int i = 0; i < lnt; i++) {
    distParsi = distPars[i];
    arma::vec extrapolParsi = extrapolPars[i];
    pbs[i] = dCount_allProbs_scalar_user(x[i], distParsi, extrapolParsi, survR,
					 nsteps, time, extrap, logFlag);
  }

  return(pbs);
}
