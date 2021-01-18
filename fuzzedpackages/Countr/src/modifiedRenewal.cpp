#include <iostream>
#include "RcppArmadillo.h"
#include "../inst/include/Countr_types.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec getRenewalExtrapolPars(arma::vec extrapolPars,
				 arma::vec extrapolPars0) {
  
  double x11 = extrapolPars[0];
  double x12 = extrapolPars[1];
  
  double x01 = extrapolPars0[0];
  double x02 = extrapolPars0[1];

  arma::vec res(2);
  res(0) = min(x11, x01);
  double bb = max(x11, x01);
  double min2 = min(x12, x02);
  res(1) = min(bb, min2);

  return(res);
}

arma::vec dePril(unsigned xnum, const arma::vec& p0, arma::vec& q,
		 arma::vec& q0, const unsigned& nsteps) {
  
  // allocate memory
  unsigned lnt = q.n_elem;
  arma::vec dfn(lnt, fill::zeros);
  arma::vec probs(2, fill::zeros);
  double nummin = (double) xnum - 1.0;
  double xmp = (double) xnum;
  unsigned n, j;
  double dn, dj, temp, fqj, s1, s2, q0n;

  if (xnum > 1) {
    dfn(0) = pow(q(0), nummin);
    for(n = 1; n <= nsteps; n++) {
      dn = (double) n;
      q0n = q(0) * dn;
      s1 = 0.0;
      s2 = 0.0;
      for(j = 1; j <= n; j++) {
	dj = (double) j;
	fqj = q(j) * dfn(n - j);
	s1 += fqj * dj;
	s2 += fqj * dn;
      }
      dfn(n) += (xmp * s1 - s2) / q0n; 
    }
   
    
    // final convolution with q0
    for(unsigned nn = 0; nn <= nsteps; nn ++) {
      n = nsteps - nn;
      temp = 0.0;
      for(j = 0; j <= n; j++) {
	temp += dfn(j) * q0(n - j);
      }
      dfn(n) = temp;
    }
  } else 
    dfn = q0;

  
  for(n = 0;n < nsteps; n++) {
    probs(0) = probs(0) + dfn(n) * p0(nsteps - n); 
    probs(1) = probs(1) + dfn(n);
  }

  bool isEven = (xnum % 2) == 0;
  if (isEven) {
    probs(0) = probs(0) + 0.5 * dfn(nsteps); 
    probs(1) = probs(1) + 0.5 * dfn(nsteps);
  }
  return(probs);
}

arma::vec doOneConvolution_dePril_odd(unsigned xnum, arma::vec& p0, 
				      arma::vec& pdf, arma::vec& pdf0, 
				      const arma::vec& fwork,
				      const arma::vec& fwork0,
				      const unsigned& nsteps, 
				      const unsigned fact) {
  double stl = 1.0;
  double stl0 = 1.0;
  double sth, sth0;
  unsigned ik;
  double xnum2 = xnum / 2;
  unsigned fact2 = fact * 0.5;

    for(unsigned i = 1 ;i <= nsteps; i++) {
      ik = fact * i;
      sth = fwork(ik);
      sth0 = fwork0(ik);
      pdf(i - 1) = stl - sth;
      pdf0(i - 1) = stl0 - sth0;
      stl = sth;
      stl0 = sth0;
      p0(i) = fwork(ik - fact2);      
    }
    
    return(dePril(xnum, p0, pdf, pdf0, nsteps - xnum2));
}

arma::vec doOneConvolution_dePril_even(unsigned xnum, arma::vec& p0, 
				      arma::vec& pdf, arma::vec& pdf0, 
				      const arma::vec& fwork,
				      const arma::vec& fwork0,
				      const unsigned& nsteps, 
				      const unsigned fact) {
  double stl = 1.0;
  double sth, sth0;
  double stl0 = 1.0;
  unsigned ik;
  double xnum2 = xnum / 2;

    for(unsigned i = 1 ;i <= nsteps; i++) {
      ik = fact * i;
      sth = fwork(ik);
      sth0 = fwork0(ik);
      pdf(i - 1) = stl - sth;
      pdf0(i - 1) = stl0 - sth0;
      stl = sth;
      stl0 = sth0;
      p0(i) = sth;      
    }

    return(dePril(xnum, p0, pdf, pdf0, nsteps - xnum2));
}

arma::vec getProbs_weibull_dePril_even(unsigned xnum,
				       const Rcpp::List distPars,
				       const std::string dist,
				       const Rcpp::List distPars0,
				       const std::string dist0,
				       arma::vec extrapolPars,
				       const unsigned& nsteps = 100, 
				       double time = 1.0, bool extrap = true) {

  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);
  
  double stl = 1.0;
  double stl0 = 1.0;
  double en = (double) nsteps;
  double xi, th, sth, sth0, tee;
  unsigned i; 
  vec pdf, p0, pdf0;
  
  if (xnum == 0) {
    probs(0) = surv(time, distPars0, dist0);
    probs(1) = 1.0;
    return(probs);
  }

  if (extrap) {
    // unsigned i8, i4, i2;
    // define the steps needed
    unsigned nsteps1 = nsteps / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    // allocate memory to the different arrays
    p0 = zeros<vec> (needed + 1);
    pdf = zeros<vec> (needed + 1);
    pdf0 = zeros<vec> (needed + 1);
    vec fwork0(needed + 1, fill::zeros);
    vec fwork(needed + 1, fill::zeros);
    
    for (i = 2; i<= needed; i += 2) {
      xi = double (i);
      tee = time * xi / en;
      fwork(i) = surv(tee, distPars, dist);;
      fwork0(i) = surv(tee, distPars0, dist0);;
    }
    
    // =========================== concolutions ================================
    // ---------- conv1 
    vec probs1 = doOneConvolution_dePril_even(xnum, p0, pdf, pdf0, 
					      fwork, fwork0, nsteps1, 8); 
    // ---------- conv2 
    vec probs2 = doOneConvolution_dePril_even(xnum, p0, pdf, pdf0, 
					      fwork, fwork0, nsteps2, 4); 
    // ---------- conv3 
    probs = doOneConvolution_dePril_even(xnum, p0, pdf, pdf0, 
					 fwork, fwork0, nsteps3, 2); 
    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4 
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
  } else {
    // just do one convolution
    p0 = zeros<vec> (nsteps + 1);
    pdf = zeros<vec> (nsteps + 1);
    pdf0 = zeros<vec> (nsteps + 1);
    
    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      th = time * xi / en;
      sth = surv(th, distPars, dist);
      sth0 = surv(th, distPars0, dist0);
      pdf(i - 1) = stl - sth;
      pdf0(i - 1) = stl0 - sth0;
      p0(i) = sth;
      stl = sth;
      stl0 = sth0;
    }
    
    probs = dePril(xnum, p0, pdf, pdf0, nsteps - xnum / 2);
  }
  return(probs);
}

arma::vec getProbs_weibull_dePril_even(unsigned xnum,
				       const Rcpp::List distPars,
				       Rcpp::Function survR,
				       const Rcpp::List distPars0,
				       Rcpp::Function survR0,
				       arma::vec extrapolPars,
				       const unsigned& nsteps = 100, 
				       double time = 1.0, bool extrap = true) {

  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);
  
  double stl = 1.0;
  double stl0 = 1.0;
  double en = (double) nsteps;
  double xi, th, sth, sth0, tee;
  unsigned i; 
  vec pdf, p0, pdf0;
  Rcpp::NumericVector rTemp;
  
  if (xnum == 0) {
    rTemp = survR0(time, distPars0);
    probs(0) = rTemp[0];;
    probs(1) = 1.0;
    return(probs);
  }

  if (extrap) {
    // unsigned i8, i4, i2;
    // define the steps needed
    unsigned nsteps1 = nsteps / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    // allocate memory to the different arrays
    p0 = zeros<vec> (needed + 1);
    pdf = zeros<vec> (needed + 1);
    pdf0 = zeros<vec> (needed + 1);
    vec fwork0(needed + 1, fill::zeros);
    vec fwork(needed + 1, fill::zeros);
    
    for (i = 2; i<= needed; i += 2) {
      xi = double (i);
      tee = time * xi / en;
      rTemp = survR(tee, distPars);
      fwork(i) = rTemp[0];
      rTemp = survR0(tee, distPars0);
      fwork0(i) = rTemp[0];
    }
    
    // =========================== concolutions ================================
    // ---------- conv1 
    vec probs1 = doOneConvolution_dePril_even(xnum, p0, pdf, pdf0, 
					      fwork, fwork0, nsteps1, 8); 
    // ---------- conv2 
    vec probs2 = doOneConvolution_dePril_even(xnum, p0, pdf, pdf0, 
					      fwork, fwork0, nsteps2, 4); 
    // ---------- conv3 
    probs = doOneConvolution_dePril_even(xnum, p0, pdf, pdf0, 
					 fwork, fwork0, nsteps3, 2); 
    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4 
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
  } else {
    // just do one convolution
    p0 = zeros<vec> (nsteps + 1);
    pdf = zeros<vec> (nsteps + 1);
    pdf0 = zeros<vec> (nsteps + 1);
    
    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      th = time * xi / en;
      rTemp = survR(th, distPars);
      sth = rTemp[0];
      rTemp = survR(th, distPars0);
      sth0 = rTemp[0];
      pdf(i - 1) = stl - sth;
      pdf0(i - 1) = stl0 - sth0;
      p0(i) = sth;
      stl = sth;
      stl0 = sth0;
    }
    
    probs = dePril(xnum, p0, pdf, pdf0, nsteps - xnum / 2);
  }
  return(probs);
}

arma::vec getProbs_weibull_dePril_odd(unsigned xnum,
				      const Rcpp::List distPars,
				      const std::string dist,
				      const Rcpp::List distPars0,
				      const std::string dist0,
				      arma::vec extrapolPars,
				      const unsigned& nsteps = 100, 
				      double time = 1.0, bool extrap = true) {
  
  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);
  
  double stl = 1.0;
  double stl0 = 1.0;
  double en = (double) nsteps;
  double xi, th, sth, sth0, tee;
  unsigned i; 
  vec pdf, p0, pdf0;
  
  if (extrap) {
    // unsigned i8, i4, i2;
    // define the steps needed
    unsigned nsteps1 = nsteps / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    // allocate memory to the different arrays
    p0 = zeros<vec> (needed + 1);
    pdf = zeros<vec> (needed + 1);
    pdf0 = zeros<vec> (needed + 1);
    vec fwork0(needed + 1, fill::zeros);
    vec fwork(needed + 1, fill::zeros);
    
    for(i = 1; i <= needed; i++) {
      xi = double (i);
      tee = time * xi / en;
      fwork(i) = surv(tee, distPars, dist);
      fwork0(i) = surv(tee, distPars0, dist0);
    }

    // =========================== concolutions ================================
    // ---------- conv1 
    // ---------- conv1 
    vec probs1 = doOneConvolution_dePril_odd(xnum, p0, pdf, pdf0, 
					     fwork, fwork0, nsteps1, 8); 
    // ---------- conv2 
    vec probs2 = doOneConvolution_dePril_odd(xnum, p0, pdf, pdf0, 
					     fwork, fwork0, nsteps2, 4); 
    // ---------- conv3 
    probs = doOneConvolution_dePril_odd(xnum, p0, pdf, pdf0, 
					fwork, fwork0, nsteps3, 2);
    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4 
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
  } else {
    // just do one convolution
    p0 = zeros<vec> (nsteps + 1);
    pdf = zeros<vec> (nsteps + 1);
    pdf0 = zeros<vec> (nsteps + 1);
    
    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      th = time * xi / en;
      sth = surv(th, distPars, dist);
      sth0 = surv(th, distPars0, dist0);
      pdf(i - 1) = stl - sth;
      pdf0(i - 1) = stl0 - sth0;
      p0(i) = surv(th - 0.5 * time / en, distPars, dist);
      stl = sth;
      stl0 = sth0;
    }
    probs = dePril(xnum, p0, pdf, pdf0, nsteps - xnum /2);
  }
  return(probs);
}

arma::vec getProbs_weibull_dePril_odd(unsigned xnum,
				      const Rcpp::List distPars,
				      Rcpp::Function survR,
				      const Rcpp::List distPars0,
				      Rcpp::Function survR0,
				      arma::vec extrapolPars,
				      const unsigned& nsteps = 100, 
				      double time = 1.0, bool extrap = true) {
  
  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);
  
  double stl = 1.0;
  double stl0 = 1.0;
  double en = (double) nsteps;
  double xi, th, sth, sth0, tee;
  unsigned i; 
  vec pdf, p0, pdf0;
  Rcpp::NumericVector rTemp;
  
  if (extrap) {
    // unsigned i8, i4, i2;
    // define the steps needed
    unsigned nsteps1 = nsteps / 4;
    unsigned nsteps2 = 2 * nsteps1;
    unsigned nsteps3 = 2 * nsteps2;
    unsigned needed = 2 * nsteps3;
    // get all the survival functions @ needed
    en = (double) needed;

    // allocate memory to the different arrays
    p0 = zeros<vec> (needed + 1);
    pdf = zeros<vec> (needed + 1);
    pdf0 = zeros<vec> (needed + 1);
    vec fwork0(needed + 1, fill::zeros);
    vec fwork(needed + 1, fill::zeros);
    
    for(i = 1; i <= needed; i++) {
      xi = double (i);
      tee = time * xi / en;
      rTemp = survR(tee, distPars);
      fwork(i) = rTemp[0];
      rTemp = survR(tee, distPars0);
      fwork0(i) = rTemp[0];
    }

    // =========================== concolutions ================================
    // ---------- conv1 
    // ---------- conv1 
    vec probs1 = doOneConvolution_dePril_odd(xnum, p0, pdf, pdf0, 
					     fwork, fwork0, nsteps1, 8); 
    // ---------- conv2 
    vec probs2 = doOneConvolution_dePril_odd(xnum, p0, pdf, pdf0, 
					     fwork, fwork0, nsteps2, 4); 
    // ---------- conv3 
    probs = doOneConvolution_dePril_odd(xnum, p0, pdf, pdf0, 
					fwork, fwork0, nsteps3, 2);
    // =============== Richardson extrapolation ================================
    double xmult1 = pow(2, extrapolPars(0)); // 2^gamma1 defined in section 4 
    double xmult = pow(2,  extrapolPars(1)); // 2^2 defined in section 4
    // update the probability
    vec pt1 = (xmult * probs2 - probs1) / (xmult - 1.0);
    vec pt2 = (xmult * probs - probs2)  / (xmult - 1.0);
    probs = (xmult1 * pt2 - pt1) / (xmult1 - 1.0);
  } else {
    // just do one convolution
    p0 = zeros<vec> (nsteps + 1);
    pdf = zeros<vec> (nsteps + 1);
    pdf0 = zeros<vec> (nsteps + 1);
    
    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      th = time * xi / en;
      rTemp = survR(th, distPars);
      sth = rTemp[0];
      rTemp = survR0(th, distPars0);
      sth0 = rTemp[0];
      pdf(i - 1) = stl - sth;
      pdf0(i - 1) = stl0 - sth0;
      rTemp = survR(th - 0.5 * time / en, distPars);
      p0(i) = rTemp[0];
      stl = sth;
      stl0 = sth0;
    }
    probs = dePril(xnum, p0, pdf, pdf0, nsteps - xnum /2);
  }
  return(probs);
}

arma::vec getProbsmodified_dePril(unsigned xnum,
				  const Rcpp::List distPars,
				  Rcpp::Function survR,
				  const Rcpp::List distPars0,
				  Rcpp::Function survR0,
				  arma::vec extrapolPars,
				  const unsigned& nsteps = 100,
				  double time = 1.0, bool extrap = true) {
  
  if ((xnum % 2) == 0)
    return(getProbs_weibull_dePril_even(xnum, distPars, survR, distPars0,
					survR0, extrapolPars,
					nsteps, time, extrap));
  else
    return(getProbs_weibull_dePril_odd(xnum, distPars, survR, distPars0,
				       survR0, extrapolPars,
				       nsteps, time, extrap));
}

// [[Rcpp::export]]
arma::vec getProbsmodified_dePril(unsigned xnum,
				  const Rcpp::List distPars,
				  const std::string dist,
				  const Rcpp::List distPars0,
				  const std::string dist0,
				  arma::vec extrapolPars,
				  const unsigned& nsteps = 100,
				  double time = 1.0, bool extrap = true) {
  
  if ((xnum % 2) == 0)
    return(getProbs_weibull_dePril_even(xnum, distPars, dist, distPars0,
					dist0, extrapolPars,
					nsteps, time, extrap));
  else
    return(getProbs_weibull_dePril_odd(xnum, distPars, dist, distPars0,
				       dist0, extrapolPars,
				       nsteps, time, extrap));
}

//' Compute count probabilities based on modified renewal process (bi)
//'
//' Compute count probabilities based on modified renewal process using
//' dePril algorithm.
//' \code{dmodifiedCount_bi} does it for the builtin distributions.
//'
//' For the modified renewal process the first arrival is allowed to have
//' a different distribution from the  time between subsequent arrivals.
//' The renewal assumption is kept.
//'
//' @param x integer (vector), the desired count values.
//' @param distPars0,distPars \code{Rcpp::List} with distribution specific slots
//'     for the first arrival and the rest of the process respectively.
//' @param dist0,dist character, name of the first and following survival
//'     distributions.
//' @param cdfout TODO
//' @inheritParams dCount_allProbs_bi
//'
//' @return vector of probabilities P(x(i)) for i = 1, ..., n where n is
//'     the length of \code{x}.
//' @export
// [[Rcpp::export]]
arma::vec dmodifiedCount_bi(arma::Col<unsigned> x, const Rcpp::List distPars,
			    const std::string dist, const Rcpp::List distPars0,
			    const std::string dist0,
			    const unsigned& nsteps = 100, 
			    double time = 1.0, bool extrap = true,
			    bool cdfout = false, bool logFlag = false) {

  arma::Col<unsigned> x_unique = unique(x);
  arma::vec vals(x.n_elem, fill::zeros);
  arma::vec val;
  arma::vec extrapolPars = getextrapolPars(distPars, dist);
  arma::vec extrapolPars0 = getextrapolPars(distPars0, dist0);
  arma::vec extrapolParsNew = getRenewalExtrapolPars(extrapolPars,
						     extrapolPars0);
  
  for (unsigned k = 0; k < x_unique.n_elem; k ++) {
    uvec ind = find(x == x_unique(k));
    arma::vec valsk(ind.n_elem, fill::ones);
    val = getProbsmodified_dePril(x_unique(k), distPars, dist, distPars0, dist0,
				  extrapolParsNew, nsteps, time, extrap);
    vals.elem(ind) = valsk * val(cdfout);
  }
  
  if (logFlag)
    return(log(vals));
  else
    return(vals);
}


//' % Compute count probabilities based on modified renewal process (user)
//'
//' % Compute count probabilities based on modified renewal process using
//' % dePril algorithm.
//' \code{dmodifiedCount_user} does the same for a user specified distribution.
//' 
//' @param survR0,survR user supplied survival function; should have 
//'     signature \code{function(t, distPars)}, where \code{t} is a positive real
//'     number (the time at which the survival function is evaluated) and
//'     \code{distPars} is a list of distribution parameters. It should return a
//'     double value (first arrival and following arrivals respectively).
//' @param extrapolPars list of same length as \code{x}, where each slot is a
//'     vector of length 2 (the extrapolation values to be used) corresponding to
//'     \code{x[i]}.
//' @inheritParams dmodifiedCount_bi
//'
//' @rdname dmodifiedCount_bi
// [[Rcpp::export]]
arma::vec dmodifiedCount_user(arma::Col<unsigned> x, const Rcpp::List distPars,
			      Rcpp::Function survR, const Rcpp::List distPars0,
			      Rcpp::Function survR0, arma::vec extrapolPars,
			      const unsigned& nsteps = 100, 
			      double time = 1.0, bool extrap = true,
			      bool cdfout = false, bool logFlag = false) {

  arma::Col<unsigned> x_unique = unique(x);
  arma::vec vals(x.n_elem, fill::zeros);
  arma::vec val;
    
  for (unsigned k = 0; k < x_unique.n_elem; k ++) {
    uvec ind = find(x == x_unique(k));
    arma::vec valsk(ind.n_elem, fill::ones);
    val = getProbsmodified_dePril(x_unique(k), distPars, survR, distPars0, survR0,
				  extrapolPars, nsteps, time, extrap);
    vals.elem(ind) = valsk * val(cdfout);
  }
  
  if (logFlag)
    return(log(vals));
  else
    return(vals);
}

//' @keywords internal
// [[Rcpp::export]]
double dmodifiedCount_scalar_bi(unsigned x, const Rcpp::List distPars, 
				const std::string dist,
				const Rcpp::List distPars0,
				const std::string dist0,
				const unsigned& nsteps = 100, 
				double time = 1.0, bool extrap = true,
				bool logFlag = false) {
  
  arma::vec extrapolPars = getextrapolPars(distPars, dist);
  arma::vec extrapolPars0 = getextrapolPars(distPars0, dist0);
  arma::vec extrapolParsNew = getRenewalExtrapolPars(extrapolPars,
						     extrapolPars0);
  
  arma::vec all = getProbsmodified_dePril(x, distPars, dist,
					  distPars0, dist0,
					  extrapolParsNew,
					  nsteps, time, extrap);
  
  double out = all(0);
  if (logFlag)
    out = log(out);
  
  return(out);
}

//' @keywords internal
// [[Rcpp::export]]
double dmodifiedCount_scalar_user(unsigned x, const Rcpp::List distPars,
				  Rcpp::Function survR,
				  const Rcpp::List distPars0,
				  Rcpp::Function survR0,
				  arma::vec extrapolPars,
				  const unsigned& nsteps = 100, 
				  double time = 1.0, bool extrap = true,
				  bool cdfout = false, bool logFlag = false) {
  
  arma::vec all = getProbsmodified_dePril(x, distPars,
					  survR, distPars0, survR0,
					  extrapolPars, nsteps, time, extrap);
  
  double out = all(0);
  if (logFlag)
    out = log(out);
      
  return(out);
}
