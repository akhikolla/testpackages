#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;
using Rcpp::Rcout;

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))


/*
 *  Kw-CWG distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  alpha in [0, 1]
 *  beta >= 0
 *  gamma >= 0
 *  a >= 0
 *  b >= 0
 *
 *  z = alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta) *
 *  (
 *     (1 - exp(-(gamma*x)**beta))**(a-1) /
 *     (alpha + (1 - alpha)*exp(-(gamma*x)**beta))**(a+1)
 *  ) *
 *  (
 *    1 - (alpha**a*(1 - exp(-(gamma*x)**beta))**a) /
 *        (alpha + (1-alpha)*exp(-(gamma*x)**beta))**a
 *  )**(b-1)
 */

inline double logpdf_kwcwg(
	double x, double alpha, double beta,
	double gamma, double a, double b)
{
	// Common term in the equation
	double aux1 = exp(-pow(gamma*x,beta));

	// Here we will factor f(x) as being A * (B^(a-1)/C^(a-1)) * (1 - D/E)^(b-1)
	double A = pow(alpha,a) * beta * gamma * a * b * pow(gamma*x,beta-1) * aux1;
	double B = 1 - aux1;
	double C = alpha + (1-alpha)*aux1;
	double D = pow(alpha,a) * pow(1-aux1,a);
	double E = pow(alpha + (1-alpha)*aux1,a);

	// return A * (B**(a-1)/C**(a+1)) * (1 - D/E)**(b-1)

	return log(A) + (a-1)*log(B) - (a+1)*log(C) + (b-1)*log(1 - D/E);
}

inline double cdf_kwcwg(
	double x, double alpha, double beta,
	double gamma, double a, double b)
{
	return 1-
		pow(1-pow(alpha,a)*
			pow(
				(1-exp(-pow(gamma*x,beta))) /
				(alpha + (1-alpha)*exp(-pow(gamma*x,beta)))
			,a)
		,b);
}

inline double invcdf_kwcwg(
	double p, double alpha, double beta,
	double gamma, double a, double b)
{
	// Common term
	double aux = pow(1-pow(1-p,1/b),1/a);

	return pow(log((alpha + (1-alpha)*aux)/(alpha*(1-aux))),1/beta)/gamma;
}

// Random number generation
inline double rng_kwcwg(
	double alpha, double beta, double gamma,
	double a, double b)
{
	return invcdf_kwcwg(R::runif(0, 1), alpha, beta, gamma, a, b);
}

// [[Rcpp::export]]
NumericVector cpp_dkwcwg(
	const NumericVector& vx,
	const NumericVector& valpha,
	const NumericVector& vbeta,
	const NumericVector& vgamma,
	const NumericVector& va,
	const NumericVector& vb,
	const bool& log_prob = false
){
	
	if(std::min(
		{vx.length(), valpha.length(), vbeta.length(),
		 vgamma.length(), va.length(), vb.length()}) < 1)
	{
		return NumericVector(0);
	}

	int maxN = std::max({
		vx.length(),
		valpha.length(),
		vbeta.length(),
		vgamma.length(),
		va.length(),
		vb.length()
	});
	NumericVector p(maxN);
	
	bool throw_warning = false;

	#pragma omp parallel for
	for(int i = 0; i < maxN; i++){
		const double x = GETV(vx, i);
		const double alpha = GETV(valpha, i);
		const double beta = GETV(vbeta, i);
		const double gamma = GETV(vgamma, i);
		const double a = GETV(va, i);
		const double b = GETV(vb, i);

		#ifdef IEEE_754
		if(ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b))
			p[i] = x+alpha+beta+gamma+a+b;
		else
		#endif
		if(alpha < 0.0 || alpha > 1.0
		   || beta < 0.0
		   || gamma < 0.0
		   || a < 0.0
		   || b < 0.0)
		{
			// Concurrency will not cause problems here.
			throw_warning = true;
			p[i] = NAN;
		} else {
			p[i] = logpdf_kwcwg(x, alpha, beta, gamma, a, b);
		}
	}

	if(!log_prob)
		p = Rcpp::exp(p);
	
	if(throw_warning)
		Rcpp::warning("NaNs produced");

	return p;
}


// [[Rcpp::export]]
NumericVector cpp_pkwcwg(
	const NumericVector& vx,
	const NumericVector& valpha,
	const NumericVector& vbeta,
	const NumericVector& vgamma,
	const NumericVector& va,
	const NumericVector& vb,
	const bool& lower_tail = true,
	const bool& log_prob = false
){

	if(std::min(
		{vx.length(), valpha.length(), vbeta.length(),
		 vgamma.length(), va.length(), vb.length()}) < 1)
	{
		return NumericVector(0);
	}

	int maxN = std::max({
		vx.length(),
		valpha.length(),
		vbeta.length(),
		vgamma.length(),
		va.length(),
		vb.length()
	});
	NumericVector p(maxN);

	bool throw_warning = false;

	#pragma omp parallel for
	for (int i = 0; i < maxN; i++){
		const double x = GETV(vx, i);
		const double alpha = GETV(valpha, i);
		const double beta = GETV(vbeta, i);
		const double gamma = GETV(vgamma, i);
		const double a = GETV(va, i);
		const double b = GETV(vb, i);

		#ifdef IEEE_754
		if(ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b))
			p[i] = x+alpha+beta+gamma+a+b;
		else
		#endif
		if(alpha < 0.0 || alpha > 1.0
		   || beta < 0.0
		   || gamma < 0.0
		   || a < 0.0
		   || b < 0.0)
		{
			// Concurrency will not cause problems here.
			throw_warning = true;
			p[i] = NAN;
		} else {
			p[i] = cdf_kwcwg(x, alpha, beta, gamma, a, b);
		}
	}

	if (!lower_tail)
		p = 1.0 - p;
	
	if (log_prob)
		p = Rcpp::log(p);
	
	if (throw_warning)
		Rcpp::warning("NaNs produced");

	return p;
}


// [[Rcpp::export]]
NumericVector cpp_qkwcwg(
	const NumericVector& vp,
	const NumericVector& valpha,
	const NumericVector& vbeta,
	const NumericVector& vgamma,
	const NumericVector& va,
	const NumericVector& vb,
	const bool& lower_tail = true,
	const bool& log_prob = false
){
	if(std::min(
		{vp.length(), valpha.length(), vbeta.length(),
		 vgamma.length(), va.length(), vb.length()}) < 1)
	{
		return NumericVector(0);
	}

	int maxN = std::max({
		vp.length(),
		valpha.length(),
		vbeta.length(),
		vgamma.length(),
		va.length(),
		vb.length()
	});
	NumericVector q(maxN);
	NumericVector pp = Rcpp::clone(vp);
	
	bool throw_warning = false;

	if (log_prob)
		pp = Rcpp::exp(pp);
	
	if (!lower_tail)
		pp = 1.0 - pp;

	#pragma omp parallel for
	for (int i = 0; i < maxN; i++){
		const double p = GETV(pp, i);
		const double alpha = GETV(valpha, i);
		const double beta = GETV(vbeta, i);
		const double gamma = GETV(vgamma, i);
		const double a = GETV(va, i);
		const double b = GETV(vb, i);

		#ifdef IEEE_754
		if(ISNAN(p) || ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b))
			q[i] = p+alpha+beta+gamma+a+b;
		else
		#endif
		if(alpha < 0.0 || alpha > 1.0
		   || beta < 0.0
		   || gamma < 0.0
		   || a < 0.0
		   || b < 0.0
		   || !VALID_PROB(p))
		{
			// Concurrency wont cause problems
			throw_warning = true;
			q[i] = NAN;
		} else {
			q[i] = invcdf_kwcwg(p, alpha, beta, gamma, a, b);
		}
	}
	
	if (throw_warning)
		Rcpp::warning("NaNs produced");

	return q;
}

// [[Rcpp::export]]
NumericVector cpp_rkwcwg(
	const int& n,
	const NumericVector& valpha,
	const NumericVector& vbeta,
	const NumericVector& vgamma,
	const NumericVector& va,
	const NumericVector& vb
){
	if(std::min(
		{valpha.length(), vbeta.length(),
		 vgamma.length(), va.length(), vb.length()}) < 1)
	{
		Rcpp::warning("NAs produced");
		return NumericVector(n, NA_REAL);
	}

	NumericVector x(n);
	
	bool throw_warning = false;

	#pragma omp parallel for
	for (int i = 0; i < n; i++){
		const double alpha = GETV(valpha, i);
		const double beta = GETV(vbeta, i);
		const double gamma = GETV(vgamma, i);
		const double a = GETV(va, i);
		const double b = GETV(vb, i);

		#ifdef IEEE_754
		if(ISNAN(alpha) || ISNAN(beta) || ISNAN(gamma) || ISNAN(a) || ISNAN(b))
			x[i] = alpha+beta+gamma+a+b;
		else
		#endif
		if(alpha < 0.0 || alpha > 1.0
		   || beta < 0.0
		   || gamma < 0.0
		   || a < 0.0
		   || b < 0.0)
		{
			// Concurrency won't cause problems
			throw_warning = true;
			x[i] = NAN;
		} else {
			x[i] = rng_kwcwg(alpha, beta, gamma, a, b);
		}
	}
	
	if (throw_warning)
		Rcpp::warning("NAs produced");

	return x;
}

