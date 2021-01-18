// #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

using namespace arma;

cube prodint(const cube & dna, int nstate, int ltimes);
cube deltaNA(const cube & nev, const mat & nrisk, int nstate, int ltimes);
cube deltaNA_LY(const cube & nev, const mat & nrisk, const mat & which_compute, int nstate, int ltimes);


RcppExport SEXP gen_msm(SEXP _times,
			SEXP _entry,
			SEXP _exit,
			SEXP _from,
			SEXP _to,
			SEXP _nstate,
			SEXP _const_modif)

{

    Rcpp::NumericVector __entry(_entry), __exit(_exit), times(_times);
    Rcpp::IntegerVector __from(_from), __to(_to);
    Rcpp::IntegerMatrix __const_modif(_const_modif);

    // ProfilerStart("/tmp/gen_msm.prof");
    
    vec entry(__entry.begin(), __entry.size(), false);
    vec exit(__exit.begin(), __exit.size(), false);
    ivec from(__from.begin(), __from.size(), false);
    ivec to(__to.begin(), __to.size(), false);
    imat const_modif(__const_modif.begin(), __const_modif.nrow(), __const_modif.ncol(), false);
    
    const int lt = times.size();
    const int n = entry.size();
    const int nstate = Rcpp::as<int>(_nstate);
    
    // define the matrices we need
    mat y(lt, nstate, fill::zeros);
    
    cube nev(nstate, nstate, lt); nev.zeros();
    cube dna(nstate, nstate, lt); dna.zeros();

    for (int j = 0; j < n; ++j) {
	for (int i = 0; i < lt; ++i) {
	    if (entry(j) < times(i) && exit(j) >= times(i)) {
		y(i, from(j) - 1) += 1;
	    }
	    if (exit(j) == times(i) && to(j) != 0) {
		nev(from(j) - 1, to(j) - 1, i) += 1;
		break;
	    }
	}
    }

    // Rcpp::Rcout << "The value of y : \n" << y << "\n";

    irowvec cc = const_modif.row(0);
    if (any(cc)) {
	umat tmp = (y >= const_modif);
	mat  which_compute = conv_to<mat>::from(tmp);
	// Nelson-Aalen (the increments) Lai and Ying
	dna = deltaNA_LY(nev, y, which_compute, nstate, lt);
    }
    else {
	// Nelson-Aalen (the increments) original
	dna = deltaNA(nev, y, nstate, lt);
    }
    
    cube est = prodint(dna, nstate, lt);	
    // ProfilerStop();
    
    return Rcpp::List::create(Rcpp::Named("n.risk") = y,
			      Rcpp::Named("n.event") = nev,
			      Rcpp::Named("dna") = dna,
			      Rcpp::Named("est") = est,
			      Rcpp::Named("time") = times);

}
	

