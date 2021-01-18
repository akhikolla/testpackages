#include "RcppArmadillo.h"
#include "../inst/include/Countr_types.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]


// =============================================================================
// ------------------------- Weibull Count -------------------------------------
// =============================================================================
// 2015-12-08 new
// like alphagenOrig but avoids repeated calculations of G(cz+1)/G(z+1)

//' Matrix of alpha terms
//'
//' Matrix of alpha terms used internally by the different Weibull count
//' functions.
//'
//' It is usually advisable to compute the alpha terms a minimum number of times
//' as it may be time consuming in general. Note that the alpha terms only depend
//' on the shape (c) parameter
//'
//' @param cc numeric, shape parameter.
//' @param jrow numeric, number of rows of the alpha matrix. See formulae (11) in
//'     \emph{McShane(2008)}.
//' @param ncol numeric, number of columns of the alpha matrix. Note that the
//'     first column corresponds to \eqn{n=0}, \eqn{n} being the count value, see
//'     formulae (11) in \emph{McShane(2008)}.
//' @return \code{jrow} x \code{ncol} (lower triangular) matrix of
//'     \eqn{\alpha_j^n} terms defined in \emph{McShane(2008)}.
//' @examples
//' ## alphagen(0.994, 6, 8)
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat alphagen(double cc , unsigned jrow, unsigned ncol) {
  arma::mat alpha(jrow, ncol, fill::zeros);
  double lgam = 0.0;
  
  for (unsigned j = 0; j < jrow; j ++) {
    alpha(j, 0) = exp(lgamma(cc * j + 1) - lgam);
    lgam += log(j + 1);
  }

  for  (unsigned n = 0; n < (ncol - 1); n ++) {
    for  (unsigned j = (n + 1); j < jrow; j ++) {
      alpha(j, n + 1) = 0;
      for  (unsigned m = n; m < j; m ++) {
	alpha(j, n + 1) = alpha(j, n + 1) + alpha(m, n)
	  * alpha(j - m, 0); //was: * exp(lgamma(cc * (j - m) + 1) - lgamma(j - m + 1));
      }
    }
  }

  return(alpha);
}

template <class TYPE>
arma::vec scalarpowmatrix(double k ,  arma::Col<TYPE> M, bool ScalarBase = true) {
  arma::vec res(M.n_elem, fill::zeros);

  if (ScalarBase) {
    for  (unsigned i = 0; i < M.n_elem; i ++) 
      res(i) = pow(k, M(i));
  } else {
    for  (unsigned i = 0; i < M.n_elem; i ++) 
      res(i) = pow(M(i), k);
  }
  return(res) ;
}

//' Univariate Weibull Count Probability
//'
//' Univariate Weibull count probability computed using matrix techniques.
//'
//' \code{dWeibullCount_mat} implements formulae (11) of \emph{McShane(2008)} to
//' compute the required probabilities.  For speed, the computations are
//' implemented in C++ and of matrix computations are used whenever possible.
//' This implementation is not efficient as it recomputes the alpha
//' matrix each time, which may slow down computation (among other things).
//'
//' \code{dWeibullCount_acc} achieves a vast (several orders of magnitude) speed
//' improvement over \code{pWeibullCountOrig}. We achieve this by using Euler-van
//' Wijngaarden techniques for accelerating the convergence of alternating series
//' and tabulation of the alpha terms available in a pre-computed matrix (shipped
//' with the package).
//'
//' When computation time is an issue, we recommend the use of
//' \code{dWeibullCount_fast}. However, \code{pWeibullCountOrig} may be more
//' accurate, especially when \code{jmax} is large.
//'
//' @param scale numeric (length 1), scale parameter of the Weibull count.
//' @param shape numeric (length 1), shape parameter of the Weibull count.
//' @param x integer (vector), the desired count values.
//' @param time double, length of the observation window (defaults to 1).
//' @param logFlag logical, if TRUE, the log of the probability will be returned.
//' @param jmax integer, number of terms used to approximate the (infinite)
//'     series.
//' @return a vector of probabilities for each component of the count vector
//'     \code{x}.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullCount_mat(arma::Col<unsigned> x,
			    double shape, double scale,
			    double time = 1.0, bool logFlag = false,
			    unsigned jmax = 50) {
  
  unsigned lnt = x.n_elem;
  arma::vec tVec(lnt, fill::zeros);
  tVec.fill(time);
  arma::vec prob(lnt, fill::zeros);
  
  arma::mat alpha_all = alphagen(shape, jmax + 1, max(x) + 1);
  arma::mat alpha_data(lnt, jmax+1, fill::zeros);

  for (unsigned i = 0; i < lnt; i ++)
    alpha_data.row(i) = alpha_all.col(x(i)).t();

  arma::mat tmp(lnt, jmax+1, fill::zeros);
  arma::vec minus1(lnt, fill::zeros);
  arma::vec tVecCC = (scale * scalarpowmatrix(shape, tVec, false));
  arma::vec Term2i(lnt, fill::zeros);
  arma::vec TermF(lnt, fill::zeros);

  for (unsigned i = 0; i < jmax; i ++) {
    arma::Col <unsigned> xi = x + i;

    minus1 = scalarpowmatrix(-1, xi);
    Term2i = scalarpowmatrix(i, tVecCC, false);
    TermF = alpha_data.col(i) * exp(-lgamma(shape * i + 1));

    tmp.col(i) = minus1 % Term2i % TermF ;
  }
  
  if (logFlag)
    prob = log(sum(tmp, 1));
  else
    prob = sum(tmp, 1);

  return(prob);
}

// Univariate Weibull Count Probability (scalar)
//
// Univariate Weibull count probability computed using matrix
// techniques (scalar).
//
// @param x unisgned scalar count
// @inheritParams dWeibullCount_mat
// @rdname dWeibullCount_mat
// @return the probability of count \code{x}.
//' @keywords internal
// [[Rcpp::export]]
double dWeibullCount_mat_scalar(unsigned x,
				double shape, double scale,
				double time = 1.0, bool logFlag = false,
				unsigned jmax = 50) {

  arma::Col<unsigned> xx(1);
  xx(0) = x;
  arma::vec res = dWeibullCount_mat(xx, shape, scale, time, logFlag, jmax);
  return(res(0));
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullCount_mat_vec(arma::Col<unsigned> x,
				arma::vec shape, arma::vec scale,
				double time = 1.0, bool logFlag = false,
				unsigned jmax = 50) {
  unsigned lnt = x.n_elem;
  arma::vec pbs(lnt, fill::zeros);
  
  if (lnt != shape.n_elem)
    stop("x and shape should have same length !");

  if (lnt != scale.n_elem)
    stop("x and scale should have same length !");

  for (unsigned i = 0; i < lnt; i++) {
    pbs[i] = dWeibullCount_mat_scalar(x[i], shape[i], scale[i], time,
				      logFlag, jmax);
  }
  
  return(pbs);
  
}

// =============================================================================
// -------------------------- Weibull-Count Euler ------------------------------
// =============================================================================
// 2015-12-08 bug-fix: replacing xi by i in some places.
arma::mat alphaTerms(double scale, double shape, arma::mat alpha_all,
		     arma::Col<unsigned> x, double t = 1.0,
		     unsigned jmax = 50) {
  double ltc = scale * pow(t, shape);
  double coeff;
  unsigned lnt = x.n_elem;
  arma::mat terms(jmax, lnt, fill::zeros);

  if(max(x) >= alpha_all.n_cols)
    stop("alpha_all does not contain enough columns!");
  
  if(jmax + max(x) > alpha_all.n_rows)
     stop("alpha_all does not contain enough rows!");

  
  for (unsigned i = 0; i < lnt; i++) {
    unsigned xi = x(i);
    arma::vec alpha = alpha_all.col(xi); // assumes xi < ncol alpha_all

    // simplified by Georgi: in McShane, eq. (11), the first term is for j=n for
    // which (-1)^(j+n) = (-1)^(2n) = 1. So, we always start with coef = 1.
    coeff = 1.0;
    for (unsigned j = xi; j < xi + jmax; j ++) {
      terms(j - xi, i)  = coeff * pow(ltc, j) * alpha(j) *
	exp(-lgamma(shape * j + 1));
      coeff = -coeff;
    }
  }

  return(terms);
}

/*
Convergence acceleration of an alternating series by the Euler transformation.
Initialize by calling the constructor with arguments nmax, an upper bound on
the number of terms to be summed, and eps, the desired accuracy.
Then make successive calls to the function next
*/
struct Eulsum {
  arma::vec wksp;
  // int n, ncv;
  unsigned n;
  int ncv;
  bool cnvgd;
  double sum,eps,lastval,lasteps;
  Eulsum(int nmax, double epss) : wksp(nmax), n(0), ncv(0),
				  cnvgd(0), sum(0.), eps(epss), lastval(0.) {}

  double next(const double term) {
    unsigned j;
    double tmp,dum;
    if (n + 1 > wksp.size())
      throw("wksp too small in eulsum");

    if (n == 0) { //Initialize:
      sum = 0.5 * (wksp[n++] = term); //Return first estimate.
    } else {
      tmp = wksp[0];
      wksp[0] = term; //Update saved quantities by van Wijngaarden’s algorithm.
      for (j=1; j < n; j++) {
	dum = wksp[j];
	wksp[j] = 0.5 * (wksp[j - 1] + tmp);
	tmp = dum;
      }
      //Favorable to increase p, and the table becomes longer.
      wksp[n] = 0.5 * (wksp[n - 1] + tmp);
      if (abs(wksp[n]) <= abs(wksp[n - 1]))
	sum += (0.5 * wksp[n++]);
      else
	//Favorable to increase n, the table doesn’t become longer.
	sum += wksp[n];
    }

    lasteps = abs(sum - lastval);
    if (lasteps <= eps) ncv++;
    if (ncv >= 2) cnvgd = 1;
    return (lastval = sum);
  }
};


//' Fast Univariate Weibull Count Probability
//'
//' @param nmax integer, an upper bound on the number of terms to be summed in
//'     the Euler-van Wijngaarden sum; default is 300 terms.
//' @param eps numeric, the desired accuracy to declare convergence.
//' @param printa logical, if \code{TRUE} print information about convergence.
//' @rdname dWeibullCount_mat
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullCount_acc(arma::Col<unsigned> x, double shape, double scale, 
			    double time = 1.0,
			    bool logFlag = false, unsigned jmax = 50,
			    int nmax = 300, double eps = 1e-10,
			    bool printa = false) {
  
  arma::mat alpha_all = alphagen(shape, jmax + max(x) + 1, max(x) + 1);
  
  arma::Col<unsigned> x_unique = unique(x);
  arma::mat terms = alphaTerms(scale, shape, alpha_all, x_unique, time, jmax);
  arma::vec vals(x.n_elem, fill::zeros);
  double val = 0.0;
  
  for (unsigned k = 0; k < x_unique.n_elem; k ++) {
    Eulsum eulsum = Eulsum(nmax, eps);
    arma::vec termsk = terms.col(k);
    unsigned i = 0;
    while (!eulsum.cnvgd && i < termsk.n_elem) {
      val = eulsum.next(termsk(i));
      i += 1;
    }
    
    if (printa) {
      if (!eulsum.cnvgd)
	Rprintf("sum did not converge !");
      else
	Rprintf(" iterations were used to reach convergence !");
    }

    uvec ind = find(x == x_unique(k));
    arma::vec valsk(ind.n_elem, fill::ones);
    vals.elem(ind) = valsk * val;
  }

  if (logFlag)
    return(log(vals));
  else
    return(vals);
}

double dWeibullCount_acc_scalar(unsigned x, double shape, double scale, 
				double time = 1.0,
				bool logFlag = false, unsigned jmax = 50,
				int nmax = 300, double eps = 1e-10,
				bool printa = false) {
  
  arma::mat alpha_all = alphagen(shape, jmax + x + 1, x + 1);
  arma::Col<unsigned> xx(1);
  xx(0) = x;
  arma::mat terms = alphaTerms(scale, shape, alpha_all, xx, time, jmax);
  double val = 0.0;
  
  Eulsum eulsum = Eulsum(nmax, eps);
  arma::vec termsk = terms.col(0);
  unsigned i = 0;
  while (!eulsum.cnvgd && i < termsk.n_elem) {
    val = eulsum.next(termsk(i));
    i += 1;
  }
    
  if (printa) {
    if (!eulsum.cnvgd)
      Rprintf("sum did not converge !");
    else
      Rprintf(" iterations were used to reach convergence !");
  }
  
  if (logFlag)
    return(log(val));
  else
    return(val);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullCount_acc_vec(arma::Col<unsigned> x,
				arma::vec shape, arma::vec scale, 
				double time = 1.0,
				bool logFlag = false, unsigned jmax = 50,
				int nmax = 300, double eps = 1e-10,
				bool printa = false) {

  unsigned lnt = x.n_elem;
  arma::vec pbs(lnt, fill::zeros);
  
  if (lnt != shape.n_elem)
    stop("x and shape should have same length !");

  if (lnt != scale.n_elem)
    stop("x and scale should have same length !");

  for (unsigned i = 0; i < lnt; i++) {
    pbs[i] = dWeibullCount_acc_scalar(x[i], shape[i], scale[i], time,
				      logFlag, jmax, nmax, eps, printa);
  }
  
  return(pbs); 

}

/*
================================================================================
---------------------------- local fast functions ------------------------------
================================================================================
*/

arma::mat alphaTerms(arma::vec scale, arma::vec shape,
		     arma::mat alpha_all, arma::Col<unsigned> x,
		     double t = 1.0, unsigned jmax = 100) {
  double coeff;
  unsigned lnt = x.n_elem;
  arma::mat terms(jmax, lnt, fill::zeros);

  if(max(x) >= alpha_all.n_cols)
    stop("alpha_all does not contain enough columns!");
  
  if(jmax + max(x) > alpha_all.n_rows)
     stop("alpha_all does not contain enough rows!");

  
  for (unsigned i = 0; i < lnt; i++) {
    unsigned xi = x(i);
    arma::vec alpha = alpha_all.col(xi); // assumes xi < ncol alpha_all
    double ltci = pow(t, shape(i)) * scale(i);
    // simplified by Georgi: in McShane, eq. (11), the first term is for j=n for
    // which (-1)^(j+n) = (-1)^(2n) = 1. So, we always start with coef = 1.
    coeff = 1.0;
    for (unsigned j = xi; j < xi + jmax; j ++) {
      terms(j - xi, i)  = coeff * pow(ltci, j) * alpha(j) *
	exp(-lgamma(shape(i) * j + 1));
      coeff = -coeff;
    }
  }

  return(terms);
}

// shape: double here
arma::mat alphaTerms(arma::vec scale, double shape,
		     arma::mat alpha_all, arma::Col<unsigned> x,
		     double t = 1.0, unsigned jmax = 100) {
  double ltc = pow(t, shape);
  double coeff;
  unsigned lnt = x.n_elem;
  arma::mat terms(jmax, lnt, fill::zeros);

  if(max(x) >= alpha_all.n_cols)
    stop("alpha_all does not contain enough columns!");
  
  if(jmax + max(x) > alpha_all.n_rows)
     stop("alpha_all does not contain enough rows!");

  
  for (unsigned i = 0; i < lnt; i++) {
    unsigned xi = x(i);
    arma::vec alpha = alpha_all.col(xi); // assumes xi < ncol alpha_all
    double ltci = ltc * scale(i);
    // simplified by Georgi: in McShane, eq. (11), the first term is for j=n for
    // which (-1)^(j+n) = (-1)^(2n) = 1. So, we always start with coef = 1.
    coeff = 1.0;
    for (unsigned j = xi; j < xi + jmax; j ++) {
      terms(j - xi, i)  = coeff * pow(ltci, j) * alpha(j) *
	exp(-lgamma(shape * j + 1));
      coeff = -coeff;
    }
  }

  return(terms);
}

// local function used in the vectorial version
// shape and scale both vectors
arma::vec dWeibullCount_fast0(arma::Col<unsigned> x, arma::vec shape, 
			      arma::vec scale, arma::mat alpha_all,
			      double t = 1.0,
			      bool logFlag = false,
			      unsigned jmax = 50, int nmax = 300,
			      double eps = 1e-10,
			      bool printa = false) {
  
  arma::mat terms = alphaTerms(scale, shape, alpha_all, x, t, jmax);
  arma::vec vals(x.n_elem, fill::zeros);
  double val = 0.0;
  
  for (unsigned k = 0; k < x.n_elem; k ++) {
    Eulsum eulsum = Eulsum(nmax, eps);
    arma::vec termsk = terms.col(k);
    unsigned i = 0;
    while (!eulsum.cnvgd && i < termsk.n_elem) {
      val = eulsum.next(termsk(i));
      i += 1;
    }
    
    if (printa) {
      if (!eulsum.cnvgd)
	Rprintf("sum did not converge !");
      else
	Rprintf(" iterations were used to reach convergence !");
    }

    vals(k) =  val;
    
  }
  
  if (logFlag)
    return(log(vals));
  else
    return(vals);
}

// local function used in the vectorial version
// double shape, vector scale
arma::vec dWeibullCount_fast0(arma::Col<unsigned> x, double shape, 
			      arma::vec scale, arma::mat alpha_all,
			      double t = 1.0,
			      bool logFlag = false,
			      unsigned jmax = 50, int nmax = 300,
			      double eps = 1e-10,
			      bool printa = false) {
  
  arma::mat terms = alphaTerms(scale, shape, alpha_all, x, t, jmax);
  arma::vec vals(x.n_elem, fill::zeros);
  double val = 0.0;
  
  for (unsigned k = 0; k < x.n_elem; k ++) {
    Eulsum eulsum = Eulsum(nmax, eps);
    arma::vec termsk = terms.col(k);
    unsigned i = 0;
    while (!eulsum.cnvgd && i < termsk.n_elem) {
      val = eulsum.next(termsk(i));
      i += 1;
    }
    
    if (printa) {
      if (!eulsum.cnvgd)
	Rprintf("sum did not converge !");
      else
	Rprintf(" iterations were used to reach convergence !");
    }

    vals(k) =  val;
    
  }
  
  if (logFlag)
    return(log(vals));
  else
    return(vals);
}

// double (scale) instead of vector scale
// shape ans scale both double
arma::vec dWeibullCount_fast0(arma::Col<unsigned> x, double shape, double scale, 
			      arma::mat alpha_all, double t = 1.0,
			      bool logFlag = false,
			      unsigned jmax = 50, int nmax = 300,
			      double eps = 1e-10,
			      bool printa = false) {
  
  arma::mat terms = alphaTerms(scale, shape, alpha_all, x, t, jmax);
  arma::vec vals(x.n_elem, fill::zeros);
  double val = 0.0;

  if (shape == 1.0) {
    for (unsigned k = 0; k < x.n_elem; k ++)
      vals(k) = R::dpois(x(k), scale, logFlag);

    return(vals);
  } else {
    for (unsigned k = 0; k < x.n_elem; k ++) {
      Eulsum eulsum = Eulsum(nmax, eps);
      arma::vec termsk = terms.col(k);
      unsigned i = 0;
      while (!eulsum.cnvgd && i < termsk.n_elem) {
	val = eulsum.next(termsk(i));
	i += 1;
      }
      
      if (printa) {
	if (!eulsum.cnvgd)
	  Rprintf("sum did not converge !");
	else
	  Rprintf(" iterations were used to reach convergence !");
      }
      
      vals(k) =  val; 
    }
    
    if (logFlag)
      return(log(vals));
    else
      return(vals);
  }
}

arma::vec dWeibullCount_fast0(unsigned x, double shape, double scale, 
			      arma::mat alpha_all,
			      double t = 1.0,
			      bool logFlag = false,
			      unsigned jmax = 50, int nmax = 300,
			      double eps = 1e-10,
			      bool printa = false) {

  arma::Col<unsigned> xx(1);
  xx(0) = x;
  return(dWeibullCount_fast0(xx, shape, scale, alpha_all, t,
			     logFlag, jmax, nmax, eps, printa)
	 );
}

// =============================================================================
// --------------------- Bivariate-Count-series specific -----------------------
// =============================================================================

// a double vector of length 2
// res(0) = cdf(y-1)
// res(1) = cdf(y)
arma::vec cdfWeibullCount(unsigned y, double shape, double scale, 
			  arma::mat alpha_all, double t = 1.0, unsigned jmax = 50,
			  int nmax = 300, double eps = 1e-10) {
  arma::vec res(2, fill::zeros);
  if (y == 0) {
    arma::Col <unsigned> x (1, fill::zeros) ;
    arma::vec proba = dWeibullCount_fast0(x, shape, scale, alpha_all, t, false, jmax,
					  nmax, eps);
    res(1) = proba(0); // res(0) is already zero
  } else {
    arma::Col <unsigned> x = seq0(y);
    arma::vec proba = dWeibullCount_fast0(x, shape, scale, alpha_all, t, false, jmax,
					  nmax, eps);
    res(0) = sum(proba(span(0, y - 1)));
    res(1) = sum(proba);
  }

  return(res);
}
