//The code in this file is taken from RcppDE, version 0.1.0, which is GPL >= 2

#ifndef _Rmalschains_evaluate_h_
#define _Rmalschains_evaluate_h_

#include <Rcpp.h>

	class EvalBase {
	public:
	    EvalBase() : neval(0) {};
	    virtual double eval(SEXP par) = 0;
	    unsigned long getNbEvals() { return neval; }
        protected:
            unsigned long int neval;
	};

	class EvalStandard : public EvalBase {
	public:
	    EvalStandard(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {} 
	    double eval(SEXP par) {
		neval++;
		return defaultfun(par);
	    }
	private:
	    SEXP fcall, env;
	    double defaultfun(SEXP par) { 		// essentialy the same as the old evaluate
		SEXP fn = PROTECT(::Rf_lang2(fcall, par)); 	// this could be done with Rcpp 
		SEXP sexp_fvec = ::Rf_eval(fn, env);	// but is still a lot slower right now
		UNPROTECT(1);		

		double f_result = REAL(sexp_fvec)[0];
		if (ISNAN(f_result)) 
		    ::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
		return(f_result); 
	    }
	};

	typedef double (*funcPtr)(SEXP, SEXP);
	class EvalCompiled : public EvalBase {
	public:
		EvalCompiled( Rcpp::XPtr<funcPtr> xptr, SEXP env_ ) {
		funptr = *(xptr);
		env = env_;
	    };
		EvalCompiled( SEXP xps, SEXP env_ ) {
		Rcpp::XPtr<funcPtr> xptr(xps);
		funptr = *(xptr);
		env = env_;
	    };
	    double eval(SEXP par) {
		neval++;
		return funptr(par, env);
	    }
	private:
	    funcPtr funptr;
	    SEXP env;
	};

#endif
