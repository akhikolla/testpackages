#include <iostream>
#include "RcppArmadillo.h"
#include "../inst/include/Countr_types.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Paper refers to the The Paper by Baker & Kharrat: Discrete
// distributions from renewal processes: fast computation of probabilities ----

// =============================================================================
// ----------------------------- local functions ----------------------------
// =============================================================================
// tool function to identify the power of 2 used
// it uses left shift operators to compute efficiently powers of 2
void whichPower(int xnum, long& po2, int& poI) {
  while (xnum != po2) {
    po2 = (po2 << 1);
    poI ++;
  }
}

// Get binary digits of a number
//
// Get binary digits of a number stored in a vector
//
// This is an adaptation from Fortran routine \code{get_bin_digits}.
// It returns the binary decomposition of a number in a vector. The vector
// contains the indexes of the power of 2 decomposition. Note that the
// Fortran table starts at index 1, when in C++ indexes starts at 0. We made
// the choice to start the result table at 0 here and adapt the rest of
// the code accordingly.
// @param xnum integer length 1 to be decomposed in power of 2
// @return an integer vector. The vector has as many elements as the
// number of power 2 indexes in the decomposition. For example, if xnum =5,
// the result will be 0, 2 (i.e 5 = 2^0 + 2^2) and if xnum = 21, the result
// will be 0, 2, 4 (i.e 21 = 2^0 + 2^2 + 2^4)
// @examples
// get_bin_digits(3) # c(0, 1)
// get_bin_digits(5) # c(0, 2)
// get_bin_digits(21) # c(0, 2, 4)
IntegerVector get_bin_digits(int xnum) {
  std::vector <int> p2;
  long po2 = 1;
  int poI = 0;

  if(xnum == 0)
    stop("0 is not accepted in binary decomposition !");

  while ( xnum != 0 ) {
    int n = xnum & (xnum - 1);
    whichPower((xnum ^ n), po2, poI);
    p2.push_back(poI);
    xnum = n;
  }

  return wrap(p2);
}

// Get double order pdf by convolution
//
// Convolves the pdf of some number of events with itself to get he double
// order pdf.
//
// This routine convolves pdfn, the pdf of some number
// of events occurring, with itself to get the pdf of double the order.
// The symmetry here means only half the number of multiplications need be done.
// More details can be found in section 3 of the Paper.
// @param pdfn double vector (passed by reference in the c++ code)
// @param h double bin length
// @param nsteps integer number of steps used.
// @return only pdfn is updated and nothing is returned.
void doublepdf(arma::vec& pdfn, const double& h, const unsigned& nsteps) {
  unsigned k, klow, j;
  double ptemp, Tmp;

  int modk = nsteps % 2;
  for(k = nsteps; k >= 1; k --) {
    // due to symmetry, we only need of the multiplication (we will x2 later)
    klow = k/2;
    ptemp  = 0.0;

    for(j = 1; j <= klow; j++)
      ptemp = ptemp + pdfn[k - j + 1] * pdfn[j];

    ptemp = ptemp + ptemp; // correct for half multiplication
    if (modk == 1) {
      Tmp = pdfn[klow + 1];
      ptemp = ptemp + Tmp * Tmp;
    }

    pdfn[k] = ptemp * h; // don't forget to x h
    modk = 1 - modk;
  }

  //shift the probability as usual from edge to midpoint
  for(k = nsteps; k >= 1; k --)
    pdfn[k] = 0.5 * (pdfn[k] + pdfn[k - 1]);
}

// get the convolution of two pdfs of different order
//
// convolves the pdf of some number of events with another pdf to get the new
// order pdf.
//
// This routine convolves \code{p}, the pdf of some number
// of events occurring, with \code{q} (ditto) to get the resulting
// new order pdf.
// @param p double vector (passed by cte reference in the c++ code). First pdf
// @param q double vector (passed by reference in the c++ code)
// @param h double bin length
// @param nsteps integer number of steps used.
// @return only q is updated and nothing is returned.
void convtwo(const arma::vec& p, arma::vec& q, const double& h,
	     const unsigned& nsteps) {
  unsigned k, j;
  double qtemp;

  for(k = nsteps; k >= 1; k --) {
    qtemp = 0.;
    for(j = 1; j <= k; j++)
      qtemp = qtemp + p[k - j + 1] * q[j];

    q[k] = qtemp * h;
  }

  for(k = nsteps; k >= 1; k --)
    q[k] = 0.5 * (q[k] + q[k - 1]);
}

// Probability that \code{xnum} > 0 events have (or at least) have
// occured by time t
//
// Probability that \code{xnum} > 0 events have (or at least) have
// occured by time t
//
// This routine calls routines to decompose num into powers of 2,
// and so to convolve the pdf in time logarithmic in \code{xnum}.
// Finally, the pdf is convolved with the probability of 0 events
// occurring to yield the required probability. (See section 3 of the paper).
// @param p0 double vector (passed by cte reference in the c++ code).
// pdf of 0 event occuring by time t
// @param pdfn double vector (passed by reference in the c++ code) and
// updated throught this routine.
// @param h double bin length (passed by cte reference in the c++ code)
// @param msteps integer number of steps used (passed by cte reference
// in the c++ code).
// @param xnum double number of events. NOTE: 0 is not acceped here.
// @return a vector of length 2 with first element the probability of
// observing \code{xnum} events by time t and second element the probability
// of observing at least \code{xnum} events by time t. \code{pdfn} will be
// updated by this routine.
arma::vec orgconv(unsigned xnum, const arma::vec& p0, arma::vec& pdfn,
		  const unsigned& msteps, const double& h) {

  // allocate memory
  unsigned lnt = pdfn.n_elem;
  arma::vec q(lnt, fill::zeros);
  arma::vec probs(2, fill::zeros);
  arma::Col< unsigned > binarray =
    as < Col< unsigned > > (get_bin_digits(xnum));
  unsigned numdigs = binarray.n_elem;
  arma::Col< unsigned > powers(max(binarray) + 1, // check size
			       fill::zeros);
  // loop index
  unsigned k, j, idoub;
  unsigned ldoub = binarray(numdigs - 1);

  // initialize powers array
  for(k = 0; k < numdigs; k++)
    powers(binarray(k)) = 1;

  if (binarray(0) == 0) // prepare q
    q.subvec(0,  msteps) = pdfn.subvec(0,  msteps);

  if (ldoub > 0) {
    for(idoub = 1; idoub <= ldoub; idoub ++) {
      // we double the pdf at each iteratio
      doublepdf(pdfn, h, msteps); //result stored in pdfn
      if (binarray(0) == idoub) {
	// prepare q by copying pdfn
	q.subvec(0,  msteps) = pdfn.subvec(0,  msteps);
      } else if (powers(idoub) == 1) {
	// convolves pdfn with q; result stored in q
	convtwo(pdfn, q, h, msteps);
      }
    }
  }

  for(j = 1; j <= msteps; j++) {
    probs(0) = probs(0) + q(msteps - j + 1) * p0(j);
    probs(1) = probs(1) + q(msteps - j + 1);
  }

  return(probs * h);
}

arma::vec doOneConvolution(unsigned xnum, arma::vec& p0,
			   arma::vec& pdf,
			   const arma::vec& fwork, const unsigned& nsteps,
			   const unsigned fact, const double h) {
  double sth;
  unsigned i, ik;
  double stl = 1.0;
  unsigned fact2 = fact * 0.5;
  for(i = 1 ;i <= nsteps; i++) {
      ik = fact * i;
      sth = fwork(ik);
      pdf(i) = stl - sth;
      stl = sth;
      p0(i) = fwork(ik - fact2);
  }

  pdf = pdf / h;
  return(orgconv(xnum, p0, pdf, nsteps, h));
}

// Compute probability of \code{xnum} events happening by convolution
//
// Compute probability of \code{xnum} events happening by convolution
// (eventually improved by Richardson extrapolation).
//
// The routine does convolutions to produce the probability of \code{xnum}
// events happening by time \code{time} (as well as the probability of at least
// \code{xnum} events) using \code{nsteps} steps, and refines result by
// Richardson extrapolation if \code{extrap} is \code{TRUE}. See section 3
// of the papaer for more details.
// @param xnum double number of events.
// @inheritParams dCount_allProbs_bi
// @return vector of length 2 with first component being the probability of
// \code{xnum} events happening by time \code{time} and second components being
// the probability of at least \code{xnum} events happening by time \code{time}.
arma::vec getProbs(unsigned xnum, const Rcpp::List distPars,
		   arma::vec extrapolPars, const std::string dist,
		   const unsigned& nsteps = 100,
		   double time = 1.0, bool extrap = false) {

  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);

  double stl = 1.0;
  double en = (double) nsteps;
  double h = time / en;
  double xi, th, sth, tee;
  unsigned i;
  vec pdf, p0;

  if (xnum == 0) {
    probs(0) = surv(time, distPars, dist);
    probs(1) = 1.0;
    return(probs);
  }

  if (extrap) { // use Richardson extrapolation to reduce the error
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
    vec fwork(needed + 1, fill::zeros);

    for(i = 1; i <= needed; i++) {
      xi = double (i);
      tee = time * xi / en;
      fwork(i) = surv(tee, distPars, dist);
    }

    // =========================== concolutions ================================
    // ---------- conv1
    h = time / ( (double) nsteps1);
    vec probs1 = doOneConvolution(xnum, p0, pdf, fwork, nsteps1, 8, h);
    // ---------- conv2
    h = time / ( (double) nsteps2);
    vec probs2 = doOneConvolution(xnum, p0, pdf, fwork, nsteps2, 4, h);
    // ---------- conv3
    h = time / ( (double) nsteps3);
    probs = doOneConvolution(xnum, p0, pdf, fwork, nsteps3, 2, h);

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

    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      th = time * xi / en;
      sth = surv(th, distPars, dist);
      pdf(i) = stl - sth;
      // NB this half left shift is same as done before:
      // move from edge to midpoint as explained in section 2 of the paper
      p0(i) = surv(th - 0.5 * time / en, distPars, dist);
      stl = sth;
    }

    pdf = pdf / h;
    probs = orgconv(xnum, p0, pdf, nsteps, h);
  }

  return(probs);
}

arma::vec getProbs(unsigned xnum, const Rcpp::List distPars,
		   arma::vec extrapolPars, Rcpp::Function survR,
		   const unsigned& nsteps = 100,
		   double time = 1.0, bool extrap = false) {

  // allocate memory and unitialize objects
  arma::vec probs(2, fill::zeros);

  double stl = 1.0;
  double en = (double) nsteps;
  double h = time / en;
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

  if (extrap) { // use Richardson extrapolation to reduce the error
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
    vec fwork(needed + 1, fill::zeros);

    for(i = 1; i <= needed; i++) {
      xi = double (i);
      tee = time * xi / en;
      rTemp = survR(tee, distPars);
      fwork(i) = rTemp[0];
    }

    // =========================== concolutions ================================
    // ---------- conv1
    h = time / ( (double) nsteps1);
    vec probs1 = doOneConvolution(xnum, p0, pdf, fwork, nsteps1, 8, h);
    // ---------- conv2
    h = time / ( (double) nsteps2);
    vec probs2 = doOneConvolution(xnum, p0, pdf, fwork, nsteps2, 4, h);
    // ---------- conv3
    h = time / ( (double) nsteps3);
    probs = doOneConvolution(xnum, p0, pdf, fwork, nsteps3, 2, h);

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

    // prepare the starting vectors
    for(i = 1; i <= nsteps; i++) {
      xi = (double) i;
      th = time * xi / en;
      rTemp = survR(th, distPars);
      sth = rTemp[0];
      pdf(i) = stl - sth;
      // NB this half left shift is same as done before:
      // move from edge to midpoint as explained in section 2 of the paper
      rTemp = survR(th - 0.5 * time / en, distPars);
      p0(i) = rTemp[0];
      stl = sth;
    }

    pdf = pdf / h;
    probs = orgconv(xnum, p0, pdf, nsteps, h);
  }

  return(probs);
}

//' Compute count probabilities using naive convolution (bi)
//'
//' Compute count probabilities using naive convolution (section 3.1) for the
//' built-in distributions
//'
//' The routine does minimum number of convolution to compute the count
//' probability P(x) sing \code{nsteps} steps, and refines result by
//' Richardson extrapolation if \code{extrap} is \code{TRUE} using the
//' algorithm of section 3.1.
//'
//' @param cdfout logical if \code{TRUE}, the cdf will be returned instead of
//' the count probability.
//' @param inheritParams dCount_allProbs_bi
//' @return vector of probabilities P(x(i)) for i = 1, ..., n where n is
//' \code{length} of \code{x}.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_naive_bi(arma::Col<unsigned> x, const Rcpp::List distPars,
			  const std::string dist,
			  const unsigned& nsteps = 100,
			  double time = 1.0, bool extrap = true,
			  bool cdfout = false, bool logFlag = false) {

  arma::vec extrapolPars = getextrapolPars(distPars, dist);
  arma::Col<unsigned> x_unique = unique(x);
  arma::vec vals(x.n_elem, fill::zeros);
  vec val;

  for (unsigned k = 0; k < x_unique.n_elem; k ++) {
    uvec ind = find(x == x_unique(k));
    arma::vec valsk(ind.n_elem, fill::ones);
    val = getProbs(x_unique(k), distPars, extrapolPars, dist, nsteps,
		   time, extrap);
    vals.elem(ind) = valsk * val(cdfout);
  }

  if (logFlag)
    return(log(vals));
  else
    return(vals);
}

//' Compute count probabilities using naive convolution (user)
//'
//' Compute count probabilities using naive convolution (section 3.1) for the
//' user passed survival function.
//'
//' @param survR Rcpp::Function user passed survival function; should have the
//' signature \code{function(t, distPars)} where \code{t} is a real number (>0)
//' where the survival function is evaluated and \code{distPars} is a list of
//' distribution parameters. It should return a double value.
//' @inheritParams dCount_naive_bi
//' @rdname dCount_naive_bi
//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_naive_user(arma::Col<unsigned> x, const Rcpp::List distPars,
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
    val = getProbs(x_unique(k), distPars, extrapolPars, survR, nsteps,
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
double dCount_naive_scalar_bi(unsigned x, const Rcpp::List distPars,
			      const std::string dist,
			      const unsigned& nsteps = 100,
			      double time = 1.0, bool extrap = true,
			      bool logFlag = false) {

  arma::vec extrapolPars = getextrapolPars(distPars, dist);
  arma::vec all = getProbs(x, distPars, extrapolPars, dist, nsteps,
			   time, extrap);

  double out = all(0);
  if (logFlag)
    out = log(out);

  return(out);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_naive_vec_bi(arma::Col<unsigned> x, const Rcpp::List distPars,
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
    pbs[i] = dCount_naive_scalar_bi(x[i], distParsi, dist, nsteps, time,
				    extrap, logFlag);
  }

  return(pbs);
}


//' @keywords internal
// [[Rcpp::export]]
double dCount_naive_scalar_user(unsigned x, const Rcpp::List distPars,
				  arma::vec extrapolPars, Rcpp::Function survR,
				  const unsigned& nsteps = 100,
				  double time = 1.0, bool extrap = true,
				  bool logFlag = false) {

  arma::vec all = getProbs(x, distPars, extrapolPars, survR, nsteps,
			      time, extrap);

  double out = all(0);
  if (logFlag)
    out = log(out);

  return(out);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dCount_naive_vec_user(arma::Col<unsigned> x,
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
    pbs[i] = dCount_naive_scalar_user(x[i], distParsi, extrapolParsi, survR,
				      nsteps, time, extrap, logFlag);
  }

  return(pbs);
}
