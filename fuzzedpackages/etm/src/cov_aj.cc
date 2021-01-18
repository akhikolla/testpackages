#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

using namespace arma;

mat cov_dna(const mat & dna, const vec & nrisk, int d, int D);

RcppExport SEXP cov_aj(SEXP __time,
		       SEXP __est,
		       SEXP __nrisk,
		       SEXP __nevent,
		       SEXP __dna)

{

    Rcpp::NumericVector _time(__time), _est(__est), _nevent(__nevent), _dna(__dna);
    Rcpp::IntegerVector dims = _est.attr("dim");
    Rcpp::NumericMatrix _nrisk(__nrisk);
    
    const int lt = _time.size();
    const int nstate = dims(0);
    const int D = nstate * nstate; // to appease Solaris compiler

    cube est(_est.begin(), nstate, nstate, lt, false);
    cube nevent(_nevent.begin(), nstate, nstate, lt, false);
    cube dna(_dna.begin(), nstate, nstate, lt, false);
    mat nrisk(_nrisk.begin(), lt, nstate, false);
    vec time(_time.begin(), _time.size(), false);

    cube cov_etm(D, D, lt);
    cov_etm.zeros();

    mat I(nstate, nstate, fill::eye);
    mat II(D, D, fill::eye);
    mat cov_deltaNA(D, D, fill::zeros);


    // first iteration
    cov_deltaNA = cov_dna(nevent.slice(0),
			  nrisk.row(0).t(),
			  nstate,
			  D);
    
    cov_etm.slice(0) = II * cov_deltaNA * II;
    
    for (int t = 1; t < lt; ++t) {
	
	mat temp_dna(dna.slice(t).begin(), nstate, nstate, false);
	mat temp_est(est.slice(t - 1).begin(), nstate, nstate, false);
	
	cov_deltaNA = cov_dna(nevent.slice(t),
			      nrisk.row(t).t(),
			      nstate,
			      D);

	cov_etm.slice(t) = kron((I + temp_dna).t(), I) * cov_etm.slice(t - 1) * kron((I+temp_dna),I) +
	    kron(I, temp_est) * cov_deltaNA * kron(I, temp_est.t());
	
    }

    return(Rcpp::wrap(cov_etm));
}
