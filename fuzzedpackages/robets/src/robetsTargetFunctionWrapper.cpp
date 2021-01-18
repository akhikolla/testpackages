
#include <vector>
#include <string>

#include <R_ext/Error.h>

//For R's Nelder-Mead solver
#include <R_ext/Applic.h>

#include <Rcpp.h>

#include "robetsTargetFunction.h"

// This function initializes all the parameters, constructs an
// object of type RobetsTargetFunction and adds an external pointer
// to this object with name "robets.xptr"
// to the environment submitted as p_rho
//
RcppExport SEXP robetsTargetFunctionInit(SEXP p_y, SEXP p_errortype, SEXP p_trendtype,
		SEXP p_seasontype, SEXP p_damped, SEXP p_lower, SEXP p_upper,
		SEXP p_opt_crit, SEXP p_nmse, SEXP p_bounds, SEXP p_m,
		SEXP p_optAlpha, SEXP p_optBeta, SEXP p_optGamma, SEXP p_optPhi, SEXP p_optSigma0, SEXP p_optInit, SEXP p_optK,
		SEXP p_givenAlpha, SEXP p_givenBeta, SEXP p_givenGamma, SEXP p_givenPhi, SEXP p_givenSigma0, SEXP p_givenInit, SEXP p_givenK,
		SEXP p_alpha, SEXP p_beta, SEXP p_gamma, SEXP p_phi, SEXP p_sigma0, SEXP p_initstate, SEXP p_k, SEXP p_rho) {

	BEGIN_RCPP;

	RobetsTargetFunction* sp = new RobetsTargetFunction();

	std::vector<double> y = Rcpp::as< std::vector<double> >(p_y);

	int errortype = Rcpp::as<int>(p_errortype);
	int trendtype = Rcpp::as<int>(p_trendtype);
	int seasontype = Rcpp::as<int>(p_seasontype);

	bool damped = Rcpp::as<bool>(p_damped);

	std::vector<double> lower = Rcpp::as< std::vector<double> >(p_lower);
	std::vector<double> upper = Rcpp::as< std::vector<double> >(p_upper);

	std::string opt_crit = Rcpp::as<std::string>(p_opt_crit);
	int nmse = Rcpp::as<int>(p_nmse);

	std::string bounds = Rcpp::as< std::string >(p_bounds);
	int m = Rcpp::as<int>(p_m);

	bool optAlpha = Rcpp::as<bool>(p_optAlpha);
	bool optBeta = Rcpp::as<bool>(p_optBeta);
	bool optGamma = Rcpp::as<bool>(p_optGamma);
	bool optPhi = Rcpp::as<bool>(p_optPhi);

	bool givenAlpha = Rcpp::as<bool>(p_givenAlpha);
	bool givenBeta = Rcpp::as<bool>(p_givenBeta);
	bool givenGamma = Rcpp::as<bool>(p_givenGamma);
	bool givenPhi = Rcpp::as<bool>(p_givenPhi);

	double alpha = Rcpp::as<double>(p_alpha);
	double beta = Rcpp::as<double>(p_beta);
	double gamma = Rcpp::as<double>(p_gamma);
	double phi = Rcpp::as<double>(p_phi);
  
  double sigma0 = Rcpp::as<double>(p_sigma0);
  std::vector<double> initstate = Rcpp::as< std::vector<double> >(p_initstate);
  double k = Rcpp::as<double>(p_k);
  
  bool optSigma0 = Rcpp::as<bool>(p_optSigma0);
  bool optInit = Rcpp::as<bool>(p_optInit);
	bool optK = Rcpp::as<bool>(p_optK);
  
  bool givenSigma0 = Rcpp::as<bool>(p_givenSigma0);
  bool givenInit = Rcpp::as<bool>(p_givenInit);
  bool givenK = Rcpp::as<bool>(p_givenK);


	sp->init(y, errortype, trendtype, seasontype, damped, lower, upper, opt_crit,
			nmse, bounds, m, optAlpha, optBeta, optGamma, optPhi, optSigma0, optInit, optK,
			givenAlpha, givenBeta, givenGamma, givenPhi, givenSigma0, givenInit, givenK,
			alpha, beta, gamma, phi, sigma0, initstate, k);

	Rcpp::Environment e(p_rho);
	e["robets.xptr"] = Rcpp::XPtr<RobetsTargetFunction>( sp, true );

	return Rcpp::wrap(e);

	END_RCPP;
}

double targetFunctionRobetsNelderMead(int n, double *par, void *ex)
{
	RobetsTargetFunction* sp = (RobetsTargetFunction*) ex;

	sp->eval(par, n);
	return sp->getObjVal();

}


RcppExport SEXP robetsNelderMead(SEXP p_var, SEXP p_env, SEXP p_abstol,
		SEXP p_intol, SEXP p_alpha, SEXP p_beta, SEXP p_gamma,
		SEXP p_trace, SEXP p_maxit)
{

	double abstol = Rcpp::as<double>(p_abstol);
	double intol = Rcpp::as<double>(p_intol);
	double alpha = Rcpp::as<double>(p_alpha);
	double beta= Rcpp::as<double>(p_beta);
	double gamma= Rcpp::as<double>(p_gamma);
  // alpha, beta, gamma are tuning parameters of the optimization method, they are not related to exp smoothing.

	int trace = Rcpp::as<int>(p_trace);
	int maxit = Rcpp::as<int>(p_maxit);

	int fncount = 0, fail=0;
	double Fmin = 0.0;

	Rcpp::NumericVector dpar(p_var);
	Rcpp::NumericVector opar(dpar.size());

	Rcpp::Environment e(p_env);
	Rcpp::XPtr<RobetsTargetFunction> sp(e.get("robets.xptr"));

	double (*funcPtr)(int n, double *par, void *ex) = targetFunctionRobetsNelderMead;

	nmmin(dpar.size(), dpar.begin(), opar.begin(), &Fmin, funcPtr,
			&fail, abstol, intol, sp, alpha, beta, gamma, trace, &fncount, maxit);

	return Rcpp::List::create(Rcpp::Named("value") = Fmin,
			Rcpp::Named("par") = opar,
			Rcpp::Named("fail") = fail,
			Rcpp::Named("fncount") = fncount);

}



