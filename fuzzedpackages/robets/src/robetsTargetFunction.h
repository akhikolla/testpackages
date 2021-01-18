#include <vector>

#include <Rcpp.h>

extern "C" {
void cpolyroot(double *opr, double *opi, int *degree,
			double *zeror, double *zeroi, Rboolean *fail);
}

class RobetsTargetFunction {

public:

	void eval(const double* p_var, int p_var_length);
	void init(std::vector<double> & p_y, int p_errortype,
  	int p_trendtype, int p_seasontype, bool p_damped,
		std::vector<double> & p_lower, std::vector<double> & p_upper, std::string p_opt_crit,
		int p_nmse, std::string p_bounds, int p_m,
		bool p_optAlpha, bool p_optBeta, bool p_optGamma, bool p_optPhi, bool p_optSigma0 ,bool p_optInit, bool p_optK,
		bool p_givenAlpha, bool p_givenBeta, bool p_givenGamma, bool p_givenPhi, bool p_givenSigma0, bool p_givenInit, bool p_givenK,
		double alpha, double beta, double gamma, double phi, double sigma0, std::vector<double> & p_initstate, double k);
  void oneEval(std::vector<double> & p_y, int p_errortype,
    int p_trendtype, int p_seasontype, bool p_damped,
		int p_nmse, int p_m,
		double alpha, double beta, double gamma, double phi, std::vector<double> & initstate, double k);

	double getObjVal() { return(objval); };
  
  double getLik(){ return(lik);};
  double getRoblik(){ return(roblik);};
  double getTau2(){ return(tau2);};
  std::vector<double> getE(){ return(e);};
  std::vector<double> getAmse(){ return(amse);};
  std::vector<double> getState(){ return(state);};
  double computeTau2(std::vector<double>& x);


private:

  void robetscalc();
  void forecast(double& l,  double& b, double *s, double *f, int& h);
  void update(double* oldl, double* l, double* oldb, double* b, double *olds, double *s, double& y);
	bool check_params();
	bool admissible();
  double rhobiweight(double x);
  double Erho();
  double psihuber(double x);
  
  double median(std::vector<double> x);

	std::vector<double> par;
	std::vector<double> y;

	int nstate;
	int errortype;
	int trendtype;
	int seasontype;
	bool damped;
	std::vector<double> par_noopt;
	std::vector<double> lower;
	std::vector<double> upper;
	std::string opt_crit;
	int nmse;
	std::string bounds;
	int m;
	int n;

	std::vector<double> state;
	double alpha, beta, gamma, phi;
  
  double sigma0;
  std::vector<double> initstate;
  double k;

	std::vector<double> e;
	std::vector<double> amse;

	double lik, objval;
  double roblik, tau2;

	bool optAlpha, optBeta, optGamma, optPhi, givenAlpha, givenBeta, givenGamma, givenPhi;
  bool optSigma0, optInit, optK, givenSigma0, givenInit, givenK;

};
