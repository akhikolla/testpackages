#include <iostream>
#include "RcppArmadillo.h"
#include "../inst/include/Countr_types.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec dePril(unsigned xnum, const arma::vec& p0,const  arma::vec& q,
		 const arma::vec& xj, const unsigned& nsteps) {

  // allocate memory
  // unsigned lnt = q.n_elem;
  arma::vec dfn(nsteps + 1, fill::zeros);
  arma::vec probs(2, fill::zeros);
  double xmp = xnum + 1.0;
  unsigned n, j, nhalf, md, endloop;
  double s1, s2, temp, qinv, xmul, db, term;

  if (xnum == 1)
    dfn = q;
  else if (xnum == 2) {
    for (n = 0; n <= nsteps; n++) {
      temp = 0.0;
      md = (n + 1) % 2;
      nhalf = n / 2;
      if (nhalf >= md) {
	endloop = nhalf - md;
	for (j = 0; j <= endloop; j++)
	  temp += q(n - j) * q(j);
      }

      temp = temp + temp;
      if (md == 1)
	temp += q(nhalf) * q(nhalf);

      dfn(n) = temp;
    }
  } else { // xnum > 2
    qinv = 1.0 / q(0);
    xmul = xmp * qinv;
    dfn(0) = pow(q(0), xnum);
    for(n = 1; n <= nsteps; n++) {
      db = xj(n);
      s1 = 0.0;
      s2 = 0.0;
      for(j = 1; j <= n; j++) {
        term = q(j) * dfn(n - j);
	s1 += xj(j) * term;
	s2 += term;
      }
      dfn(n) = xmul * s1/db - qinv*s2;
    }
  }

  for(n = 0;n < nsteps; n++) {
    probs(0) = probs(0) + dfn(n) * p0(nsteps - n);
    probs(1) = probs(1) + dfn(n);
  }

  if ((xnum % 2) == 0)
    probs = probs + 0.5 * dfn(nsteps);

  return(probs);
}

arma::vec doOneConvolution_dePril_even(unsigned xnum, arma::vec& p0,
				       arma::vec& pdf,
				       const arma::vec& fwork,
				       const arma::vec& xj,
				       const unsigned& nsteps,
				       const unsigned fact) {
  double sth;
  unsigned i, ik;
  double stl = 1.0;
  double xnum2 = xnum / 2;
  for(i = 1 ;i <= nsteps; i++) {
      ik = fact * i;
      sth = fwork(ik);
      pdf(i - 1) = stl - sth;
      stl = sth;
      p0(i) = sth;
  }

  return(dePril(xnum, p0, pdf, xj, nsteps - xnum2));
}

arma::vec doOneConvolution_dePril_odd(unsigned xnum, arma::vec& p0,
				      arma::vec& pdf,
				      const arma::vec& fwork,
				      const arma::vec& xj,
				      const unsigned& nsteps,
				      const unsigned fact) {

  double sth;
  unsigned i, ik;
  double stl = 1.0;
  unsigned fact2 = fact * 0.5;
  for(i = 1 ;i <= nsteps; i++) {
      ik = fact * i;
      sth = fwork(ik);
      pdf(i - 1) = stl - sth;
      stl = sth;
      p0(i) = fwork(ik - fact2);
  }

  return(dePril(xnum, p0, pdf, xj, nsteps - xnum / 2));
}

arma::vec getProbs_dePril_even(unsigned xnum, const Rcpp::List distPars,
			       arma::vec extrapolPars, const std::string dist,
			       const unsigned& nsteps, double time,
			       bool extrap) {

  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);

  double stl, en;
  double xi, th, sth, tee;
  unsigned i;
  vec pdf, p0;

  if (xnum == 0) {
    probs(0) = surv(time, distPars, dist);
    probs(1) = 1.0;
    return(probs);
  }

  if (extrap) {
    // unsigned i8, i4, i2;
     unsigned adjust_steps = 0;
    // check if nsteps is enough
    if (nsteps < 2 * xnum) {
      unsigned di = 2 * xnum - nsteps;
      adjust_steps = di + 20;
    }
    
    // define the steps needed
    unsigned nsteps1 = (nsteps + adjust_steps) / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    // allocate memory to the different arrays
    p0 = zeros<vec> (needed + 1);
    pdf = zeros<vec> (needed + 1);
    vec fwork(needed + 1, fill::zeros);
    vec xjarray(needed + 1, fill::zeros);

    for (i = 2; i<= needed; i += 2) {
      xi = double (i);
      xjarray(i) = xi;
      xjarray(i - 1) = xi - 1.0;
      tee = time * xi / en;
      fwork(i) = surv(tee, distPars, dist);
    }

    // =========================== concolutions ================================
    // ---------- conv1
    vec probs1 = doOneConvolution_dePril_even(xnum, p0, pdf, fwork, xjarray,
					      nsteps1, 8);
    // ---------- conv2
    vec probs2 = doOneConvolution_dePril_even(xnum, p0, pdf, fwork, xjarray,
					      nsteps2, 4);
    // ---------- conv3
    probs = doOneConvolution_dePril_even(xnum, p0, pdf, fwork, xjarray,
					 nsteps3, 2);

    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
  } else {
    stl = 1.0;
    en = (double) nsteps;
    // just do one convolution
    p0 = zeros<vec> (nsteps + 1);
    pdf = zeros<vec> (nsteps + 1);
    vec xjarray(nsteps + 1, fill::zeros);

    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      xjarray(i) = xi;
      th = time * xi / en;
      sth = surv(th, distPars, dist);
      pdf(i - 1) = stl - sth;
      p0(i) = sth;
      stl = sth;
    }
    probs = dePril(xnum, p0, pdf, xjarray, nsteps -  xnum / 2);
  }

  return(probs);
}

arma::vec getProbs_dePril_even(unsigned xnum, const Rcpp::List distPars,
			       arma::vec extrapolPars, Rcpp::Function survR,
			       const unsigned& nsteps, double time,
			       bool extrap) {

  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);

  double stl, en;
  double xi, th, sth, tee;
  unsigned i;
  vec pdf, p0;
  Rcpp::NumericVector rTemp;

  if (xnum == 0) {
    rTemp = survR(time, distPars);
    probs(0) = rTemp[0];
    probs(1) = 1.0;
    return(probs);
  }

  if (extrap) {
    // unsigned i8, i4, i2;
     unsigned adjust_steps = 0;
    // check if nsteps is enough
    if (nsteps < 2 * xnum) {
      unsigned di = 2 * xnum - nsteps;
      adjust_steps = di + 20;
    }
    
    // define the steps needed
    unsigned nsteps1 = (nsteps + adjust_steps) / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    // allocate memory to the different arrays
    p0 = zeros<vec> (needed + 1);
    pdf = zeros<vec> (needed + 1);
    vec fwork(needed + 1, fill::zeros);
    vec xjarray(needed + 1, fill::zeros);

    for (i = 2; i<= needed; i += 2) {
      xi = double (i);
      xjarray(i) = xi;
      xjarray(i - 1) = xi - 1.0;
      tee = time * xi / en;
      rTemp = survR(tee, distPars);
      fwork(i) = rTemp[0];
    }

    // =========================== concolutions ================================
    // ---------- conv1
    vec probs1 = doOneConvolution_dePril_even(xnum, p0, pdf, fwork, xjarray,
					      nsteps1, 8);
    // ---------- conv2
    vec probs2 = doOneConvolution_dePril_even(xnum, p0, pdf, fwork, xjarray,
					      nsteps2, 4);
    // ---------- conv3
    probs = doOneConvolution_dePril_even(xnum, p0, pdf, fwork, xjarray,
					 nsteps3, 2);

    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
  } else {
    stl = 1.0;
    en = (double) nsteps;
    // just do one convolution
    p0 = zeros<vec> (nsteps + 1);
    pdf = zeros<vec> (nsteps + 1);
    vec xjarray(nsteps + 1, fill::zeros);

    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      xjarray(i) = xi;
      th = time * xi / en;
      rTemp =  survR(th, distPars);
      sth = rTemp[0];
      pdf(i - 1) = stl - sth;
      p0(i) = sth;
      stl = sth;
    }
    probs = dePril(xnum, p0, pdf, xjarray, nsteps -  xnum / 2);
  }

  return(probs);
}

arma::vec getProbs_dePril_odd(unsigned xnum, const Rcpp::List distPars,
			      arma::vec extrapolPars, const std::string dist,
			      const unsigned& nsteps, double time,
			      bool extrap) {

  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);

  double stl, en;
  double xi, th, sth, tee;
  unsigned i;
  vec pdf, p0;

  if (extrap) {
    // unsigned i8, i4, i2;
     unsigned adjust_steps = 0;
    // check if nsteps is enough
    if (nsteps < 2 * xnum) {
      unsigned di = 2 * xnum - nsteps;
      adjust_steps = di + 20;
    }
    
    // define the steps needed
    unsigned nsteps1 = (nsteps + adjust_steps) / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    // allocate memory to the different arrays
    p0 = zeros<vec> (needed + 1);
    pdf = zeros<vec> (needed + 1);
    vec fwork(needed + 1, fill::zeros);
    vec xjarray(needed + 1, fill::zeros);

    for(i = 1; i <= needed; i++) {
      xi = double (i);
      xjarray(i) = xi;
      tee = time * xi / en;
      fwork(i) =  surv(tee, distPars, dist);
    }

    // =========================== concolutions ================================
    // ---------- conv1
    vec probs1 = doOneConvolution_dePril_odd(xnum, p0, pdf, fwork, xjarray,
					     nsteps1, 8);
    // ---------- conv2
    vec probs2 = doOneConvolution_dePril_odd(xnum, p0, pdf, fwork, xjarray,
					     nsteps2, 4);
    // ---------- conv3
    probs = doOneConvolution_dePril_odd(xnum, p0, pdf, fwork, xjarray,
					nsteps3, 2);
    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
  } else {
    stl = 1.0;
    en = (double) nsteps;
    // just do one convolution
    p0 = zeros<vec> (nsteps + 1);
    pdf = zeros<vec> (nsteps + 1);
    vec xjarray(nsteps + 1, fill::zeros);

    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      xjarray(i) = xi;
      th = time * xi / en;
      sth = surv(th, distPars, dist);
      pdf(i - 1) = stl - sth;
      p0(i) = surv(th - 0.5 * time / en, distPars, dist);
      stl = sth;
    }
    probs = dePril(xnum, p0, pdf, xjarray, nsteps -  xnum / 2);
  }

  return(probs);
}

arma::vec getProbs_dePril_odd(unsigned xnum, const Rcpp::List distPars,
			      arma::vec extrapolPars, Rcpp::Function survR,
			      const unsigned& nsteps, double time,
			      bool extrap) {

  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);

  double stl, en;
  double xi, th, sth, tee;
  unsigned i;
  vec pdf, p0;
  Rcpp::NumericVector rTemp;

  if (extrap) {
    // unsigned i8, i4, i2;
     unsigned adjust_steps = 0;
    // check if nsteps is enough
    if (nsteps < 2 * xnum) {
      unsigned di = 2 * xnum - nsteps;
      adjust_steps = di + 20;
    }
    
    // define the steps needed
    unsigned nsteps1 = (nsteps + adjust_steps) / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    // allocate memory to the different arrays
    p0 = zeros<vec> (needed + 1);
    pdf = zeros<vec> (needed + 1);
    vec fwork(needed + 1, fill::zeros);
    vec xjarray(needed + 1, fill::zeros);

    for(i = 1; i <= needed; i++) {
      xi = double (i);
      xjarray(i) = xi;
      tee = time * xi / en;
      rTemp = survR(tee, distPars);
      fwork(i) =  rTemp[0];
    }

    // =========================== concolutions ================================
    // ---------- conv1
    vec probs1 = doOneConvolution_dePril_odd(xnum, p0, pdf, fwork, xjarray,
					     nsteps1, 8);
    // ---------- conv2
    vec probs2 = doOneConvolution_dePril_odd(xnum, p0, pdf, fwork, xjarray,
					     nsteps2, 4);
    // ---------- conv3
    probs = doOneConvolution_dePril_odd(xnum, p0, pdf, fwork, xjarray,
					nsteps3, 2);
    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
  } else {
    stl = 1.0;
    en = (double) nsteps;
    // just do one convolution
    p0 = zeros<vec> (nsteps + 1);
    pdf = zeros<vec> (nsteps + 1);
    vec xjarray(nsteps + 1, fill::zeros);

    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      xjarray(i) = xi;
      th = time * xi / en;
      rTemp = survR(th, distPars);
      sth = rTemp[0];
      pdf(i - 1) = stl - sth;
      rTemp = survR(th - 0.5 * time / en, distPars);
      p0(i) = rTemp[0];
      stl = sth;
    }
    probs = dePril(xnum, p0, pdf, xjarray, nsteps -  xnum / 2);
  }

  return(probs);
}

arma::vec getProbs_dePril(unsigned xnum, const Rcpp::List distPars,
			  arma::vec extrapolPars, const std::string dist,
			  const unsigned& nsteps = 100,
			  double time = 1.0, bool extrap = true) {

  if ((xnum % 2) == 0)
    return(getProbs_dePril_even(xnum, distPars, extrapolPars, dist,
					nsteps, time, extrap));
  else
    return(getProbs_dePril_odd(xnum, distPars, extrapolPars, dist,
				       nsteps, time, extrap));
}

arma::vec getProbs_dePril(unsigned xnum, const Rcpp::List distPars,
			  arma::vec extrapolPars, Rcpp::Function survR,
			  const unsigned& nsteps = 100,
			  double time = 1.0, bool extrap = true) {

  if ((xnum % 2) == 0)
    return(getProbs_dePril_even(xnum, distPars, extrapolPars, survR,
					nsteps, time, extrap));
  else
    return(getProbs_dePril_odd(xnum, distPars, extrapolPars, survR,
				       nsteps, time, extrap));
}

//' Compute count probabilities using dePril convolution (bi)
//'
//' Compute count probabilities using dePril convolution (section 3.2) for the
//' built-in distributions
//'
//' The routine does minimum number of convolution using dePril trick
//' to compute the count
//' probability P(x) sing \code{nsteps} steps, and refines result by
//' Richardson extrapolation if \code{extrap} is \code{TRUE} using the
//' algorithm of section 3.2.
//'
//' @param inheritParams dCount_naive_bi
//' @return vector of probabilities P(x(i)) for i = 1, ..., n where n is
//' \code{length} of \code{x}.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_dePril_bi(arma::Col<unsigned> x, const Rcpp::List distPars,
			   const std::string dist,
			   const unsigned& nsteps = 100,
			   double time = 1.0, bool extrap = true,
			   bool cdfout = false, bool logFlag = false) {

  arma::Col<unsigned> x_unique = unique(x);
  arma::vec vals(x.n_elem, fill::zeros);
  arma::vec val;
  arma::vec extrapolPars = getextrapolPars(distPars, dist);

  for (unsigned k = 0; k < x_unique.n_elem; k ++) {
    uvec ind = find(x == x_unique(k));
    arma::vec valsk(ind.n_elem, fill::ones);
    val = getProbs_dePril(x_unique(k), distPars, extrapolPars, dist, nsteps,
		   time, extrap);
    vals.elem(ind) = valsk * val(cdfout);
  }

  if (logFlag)
    return(log(vals));
  else
    return(vals);
}

//' Compute count probabilities using dePril convolution (user)
//'
//' Compute count probabilities using dePril convolution (section 3.2) for the
//' user passed survival function.
//'
//' @param survR Rcpp::Function user passed survival function; should have the
//' signature \code{function(t, distPars)} where \code{t} is a real number (>0)
//' where the survival function is evaluated and \code{distPars} is a list of
//' distribution parameters. It should return a double value.
//' @inheritParams dCount_dePril_bi
//' @rdname dCount_dePril_bi
//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_dePril_user(arma::Col<unsigned> x, const Rcpp::List distPars,
			    arma::vec extrapolPars, Rcpp::Function survR,
			    const unsigned& nsteps = 100,
			    double time = 1.0, bool extrap = true,
			    bool cdfout = false, bool logFlag = false) {
  arma::Col<unsigned> x_unique = unique(x);
  arma::vec vals(x.n_elem, fill::zeros);
  vec val;

  for (unsigned k = 0; k < x_unique.n_elem; k ++) {
    uvec ind = find(x == x_unique(k));
    arma::vec valsk(ind.n_elem, fill::ones);
    val = getProbs_dePril(x_unique(k), distPars, extrapolPars, survR, nsteps,
			  time, extrap);
    vals.elem(ind) = valsk * val(cdfout);
  }

  if (logFlag)
    return(log(vals));
  else
    return(vals);
}

//' @keywords internal
// [[Rcpp::export]]
double dCount_dePril_scalar_bi(unsigned x, const Rcpp::List distPars,
			       const std::string dist,
			       const unsigned& nsteps = 100,
			       double time = 1.0, bool extrap = true,
			       bool logFlag = false) {

  arma::vec extrapolPars = getextrapolPars(distPars, dist);
  arma::vec all = getProbs_dePril(x, distPars, extrapolPars, dist, nsteps,
			   time, extrap);

  double out = all(0);
  if (logFlag)
    out = log(out);

  return(out);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_dePril_vec_bi(arma::Col<unsigned> x, const Rcpp::List distPars,
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
    pbs[i] = dCount_dePril_scalar_bi(x[i], distParsi, dist, nsteps, time,
				     extrap, logFlag);
  }

  return(pbs);
}

//' @keywords internal
// [[Rcpp::export]]
double dCount_dePril_scalar_user(unsigned x, const Rcpp::List distPars,
				  arma::vec extrapolPars, Rcpp::Function survR,
				  const unsigned& nsteps = 100,
				  double time = 1.0, bool extrap = true,
				  bool logFlag = false) {

  arma::vec all = getProbs_dePril(x, distPars, extrapolPars, survR, nsteps,
			      time, extrap);

  double out = all(0);
  if (logFlag)
    out = log(out);

  return(out);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_dePril_vec_user(arma::Col<unsigned> x, const Rcpp::List distPars,
				 const Rcpp::List extrapolPars, Rcpp::Function survR,
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
    pbs[i] = dCount_dePril_scalar_user(x[i], distParsi, extrapolParsi, survR,
				       nsteps, time, extrap, logFlag);
  }

  return(pbs);
}
