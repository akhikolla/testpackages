#include "RcppArmadillo.h"
#include "../inst/include/Countr_types.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

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

//' Univariate Weibull Count Probability with gamma heterogeneity
//'
//' @param r,alpha double gamma parameters
//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullgammaCount_mat(arma::Col<unsigned> x, double shape,
				 double r, double alpha,
				 double time = 1.0, bool logFlag = false,
				 unsigned jmax = 100) {

  unsigned lnt = x.n_elem;
  arma::vec tVec(lnt, fill::zeros);
  tVec.fill(time);
  arma::vec prob(lnt, fill::zeros);
  
  arma::mat alpha_all = alphagen(shape, jmax + 1, max(x) + 1);		
  arma::mat alpha_data(lnt, jmax + 1, fill::zeros);

  for (unsigned i = 0; i < lnt; i ++)
    alpha_data.row(i) = alpha_all.col(x(i)).t();

  arma::mat tmp(lnt, jmax + 1, fill::zeros);
  arma::vec minus1(lnt, fill::zeros);
  arma::vec tVecCC = log(tVec) * shape;
  double fact = 0.0;
  
  for (unsigned j = 0; j <= jmax; j ++) {
    arma::Col <unsigned> xj = x + j;
    
    minus1 = scalarpowmatrix(-1, xj);

    tmp.col(j) = minus1 % 
      exp(
	  j * tVecCC + 
	  log(alpha_data.col(j)) -
	  lgamma(shape * j + 1)
	  +  fact);
    fact += log(r + j) - log(alpha); 
  }

  if (logFlag)
    prob = log(sum(tmp, 1));
  else
    prob = sum(tmp, 1);

  return(prob);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullgammaCount_mat_vec(arma::Col<unsigned> x,
				     arma::vec shape,
				     double r, double alpha,
				     double time = 1.0,
				     bool logFlag = false,
				     unsigned jmax = 100) {

  unsigned lnt = x.n_elem;
  arma::vec pbs(lnt, fill::zeros);
  arma::Col<unsigned> xi(1, fill::zeros);
  arma::vec pbi;
  
  if (lnt != shape.n_elem)
    stop("x and shape should have same length !");

  for (unsigned i = 0; i < lnt; i++) {
    xi[0] = x[i];
    pbi = dWeibullgammaCount_mat(xi, shape[i], r, alpha,  time,
				 logFlag, jmax);
    pbs[i] = pbi[0];
  }
  
  return(pbs);
}

//' Univariate Weibull Count Probability with gamma and
//' covariate heterogeneity
//'
//' @param r numeric shape of the gamma distribution
//' @param alpha numeric rate of the gamma distribution
//' @param Xcovar matrix covariates value
//' @param beta numeric vector of slopes
//' @param x,cc,t,logFlag,jmax TODO
//' keywords internal
// [[Rcpp::export]]
arma::vec dWeibullgammaCount_mat_Covariates(arma::Col<unsigned> x, double cc,
					    double r, double alpha,
					    arma::mat Xcovar, arma::vec beta,
					    double t = 1.0, bool logFlag = false,
					    unsigned jmax = 100) {
  unsigned lnt = x.n_elem;
  arma::vec tVec(lnt, fill::zeros);
  tVec.fill(t);
  arma::vec prob(lnt, fill::zeros);
  arma::mat alpha_all = alphagen(cc, jmax + 1, max(x) + 1);
  arma::mat alpha_data(lnt, jmax + 1, fill::zeros);

  for (unsigned i = 0; i < lnt; i ++)
    alpha_data.row(i) = alpha_all.col(x(i)).t();

  arma::mat tmp(lnt, jmax + 1, fill::zeros);
  arma::vec minus1(lnt, fill::zeros);
  arma::vec tVecCC = log(tVec)*cc +  Xcovar * beta;
  double fact = 0.0;

  for (unsigned j = 0; j <= jmax; j ++) {
    arma::Col <unsigned> xj = x + j;

    minus1 = scalarpowmatrix(-1, xj);
        
    tmp.col(j) = minus1 %
      exp(
	  j * tVecCC + 
	  log(alpha_data.col(j)) -
	  lgamma(cc * j + 1)
	  +  fact);
    fact += log(r + j) - log(alpha); 
  }

  if (logFlag)
    prob = log(sum(tmp, 1));
  else
    prob = sum(tmp, 1);

  return(prob);
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullgammaCount_mat_Covariates_vec(arma::Col<unsigned> x,
						arma::vec cc,
						double r, double alpha,
						arma::mat Xcovar, arma::vec beta,
						double t = 1.0,
						bool logFlag = false,
						unsigned jmax = 100) {
  unsigned lnt = x.n_elem;
  arma::vec pbs(lnt, fill::zeros);
  arma::Col<unsigned> xi(1, fill::zeros);
  arma::vec pbi;
  
  if (lnt != cc.n_elem)
    stop("x and shape should have same length !");

  for (unsigned i = 0; i < lnt; i++) {
    xi[0] = x[i];
    pbi = dWeibullgammaCount_mat_Covariates(xi, cc[i], r, alpha, Xcovar, beta,
					    t, logFlag, jmax);
    pbs[i] = pbi[0];
  }
  
  return(pbs);
}

// =============================================================================
// -------------------------- Weibull-Count Euler ------------------------------
// =============================================================================
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

arma::mat alphaTerms_gammaHet(double r, double alphaPar, arma::vec xb,
			      double shape,
			      arma::mat alpha_all, arma::Col<unsigned> x,
			      double t = 1.0, unsigned jmax = 100) {
  double ltc = log(t) * shape;
  double coeff, fact;
  unsigned lnt = x.n_elem;
  arma::mat terms(jmax, lnt, fill::zeros);

  if(max(x) >= alpha_all.n_cols)
    stop("alpha_all does not contain enough columns!");
  
  if(jmax + max(x) > alpha_all.n_rows)
     stop("alpha_all does not contain enough rows!");

  
  for (unsigned i = 0; i < lnt; i++) {
    unsigned xi = x(i);
    arma::vec alpha = log(alpha_all.col(xi)); // assumes xi < ncol alpha_all
    double ltci = ltc +  xb(i);
    fact = lgamma(r + xi) - lgamma(r) - log(alphaPar) * xi;
    coeff = 1.0;
    for (unsigned j = xi; j < xi + jmax; j ++) {
      terms(j - xi, i)  = coeff *
	exp(ltci *  j +
	    alpha(j) -lgamma(shape * j + 1) + fact
	    );
      coeff = -coeff;
      fact += log(r + j) - log(alphaPar);
    }
  }

  return(terms);
}

arma::mat alphaTerms_gammaHet(double r, double alphaPar, double shape,
			      arma::mat alpha_all, arma::Col<unsigned> x,
			      double t = 1.0, unsigned jmax = 100) {
  double ltc = log(t) * shape;
  double coeff, fact;
  unsigned lnt = x.n_elem;
  arma::mat terms(jmax, lnt, fill::zeros);

  if(max(x) >= alpha_all.n_cols)
    stop("alpha_all does not contain enough columns!");
  
  if(jmax + max(x) > alpha_all.n_rows)
     stop("alpha_all does not contain enough rows!");

  
  for (unsigned i = 0; i < lnt; i++) {
    unsigned xi = x(i);
    arma::vec alpha = log(alpha_all.col(xi)); // assumes xi < ncol alpha_all
    fact = lgamma(r + xi) - lgamma(r) - log(alphaPar) * xi;
    coeff = 1.0;
    for (unsigned j = xi; j < xi + jmax; j ++) {
      terms(j - xi, i)  = coeff *
	exp(j * ltc +
	    alpha(j) -lgamma(shape * j + 1) +
	    fact
	    );
      coeff = -coeff;
      fact += log(r + j) - log(alphaPar);
    }
  }

  return(terms);
}

//' Fast Univariate Weibull Count Probability with gamma heterogeneity
//'
//' @param x integer (vector), the desired count values.
//' @param shape numeric, shape parameter.
//' @param time double, length of the observation window (defaults to 1).
//' @param r numeric shape of the gamma distribution
//' @param alpha numeric rate of the gamma distribution
//' @param printa logical, if \code{TRUE} print information about convergence.
//' @param logFlag logical, if TRUE, the log of the probability will be returned.
//' @param jmax integer, number of terms used to approximate the (infinite)
//'     series.
//' @param nmax integer, an upper bound on the number of terms to be summed in
//'     the Euler-van Wijngaarden sum; default is 300 terms.
//' @param eps numeric, the desired accuracy to declare convergence.
//' @rdname dWeibullgammaCount
//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullgammaCount_acc(arma::Col<unsigned> x, double shape,
				 double r, double alpha, double time = 1.0,
				 bool logFlag = false, unsigned jmax = 100,
				 int nmax = 300, double eps = 1e-10,
				 bool printa = false) {
  
  arma::mat alpha_all = alphagen(shape, jmax + max(x) + 1, max(x) + 1);
    
  arma::Col<unsigned> x_unique = unique(x);
  arma::mat terms = alphaTerms_gammaHet(r, alpha, shape, alpha_all,
					x_unique, time, jmax);
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

//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullgammaCount_acc_vec(arma::Col<unsigned> x, arma::vec shape,
				     double r, double alpha, double time = 1.0,
				     bool logFlag = false, unsigned jmax = 100,
				     int nmax = 300, double eps = 1e-10,
				     bool printa = false) {
  unsigned lnt = x.n_elem;
  arma::vec pbs(lnt, fill::zeros);
  arma::Col<unsigned> xi(1, fill::zeros);
  arma::vec pbi;
  
  if (lnt != shape.n_elem)
    stop("x and shape should have same length !");

  for (unsigned i = 0; i < lnt; i++) {
    xi[0] = x[i];
    pbi = dWeibullgammaCount_acc(xi, shape[i], r, alpha,
				 time, logFlag, jmax, nmax, eps, printa);
    pbs[i] = pbi[0];
  }
  
  return(pbs);
}



// Fast Univariate Weibull Count Probability with gamma and covariate
// heterogeneity
//
//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullgammaCount_acc_Covariates(arma::Col<unsigned> x, double cc,
					    double r, double alpha,
					    arma::mat Xcovar, arma::vec beta,
					    double t = 1.0,
					    bool logFlag = false,
					    unsigned jmax = 100,
					    int nmax = 300, double eps = 1e-10,
					    bool printa = false) {

  if(x.n_elem != Xcovar.n_rows)
    stop("check dimesnion Xcovar and x!");
  if(beta.n_elem != Xcovar.n_cols)
    stop("check dimesnion Xcovar and beta!");
  
  arma::mat alpha_all = alphagen(cc, jmax + max(x) + 1, max(x) + 1);
  
  arma::mat terms = alphaTerms_gammaHet(r, alpha, Xcovar * beta, cc, alpha_all,
					x, t, jmax);
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

//' @keywords internal
// [[Rcpp::export]]
arma::vec dWeibullgammaCount_acc_Covariates_vec(arma::Col<unsigned> x,
						arma::vec cc,
						double r, double alpha,
						arma::mat Xcovar, arma::vec beta,
						double t = 1.0,
						bool logFlag = false,
						unsigned jmax = 100,
						int nmax = 300, double eps = 1e-10,
						bool printa = false) {
  unsigned lnt = x.n_elem;
  arma::vec pbs(lnt, fill::zeros);
  arma::Col<unsigned> xi(1, fill::zeros);
  arma::vec pbi;
  
  if (lnt != cc.n_elem)
    stop("x and cc should have same length !");

  for (unsigned i = 0; i < lnt; i++) {
    xi[0] = x[i];
    pbi = dWeibullgammaCount_acc_Covariates(xi, cc[i], r, alpha, Xcovar, beta,
					    t, logFlag, jmax, nmax, eps, printa);
    pbs[i] = pbi[0];
  }
  
  return(pbs);
}
