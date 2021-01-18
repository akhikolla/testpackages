#include <RcppArmadillo.h>

using namespace arma;

vec my_diff(const vec bouh) {

    int lbouh = bouh.size();
    vec res(lbouh - 1); res.zeros();
    
    for (int i = 0; i < (lbouh - 1); ++i) {
	res[i] = bouh[i+1] - bouh[i];
    }

    return res;
}
	    
cube prodint(const cube & dna, int nstate, int ltimes) {

    cube aj(dna.begin(), nstate, nstate, ltimes);
    mat I = eye<mat>(nstate, nstate);

    aj.slice(0) = aj.slice(0) + I;

    for (int i = 1; i < ltimes; ++i) {
	aj.slice(i) = aj.slice(i-1) * (I + aj.slice(i));
    }

    return aj;
}


cube deltaNA(const cube & nev, const mat & nrisk, int nstate, int ltimes) {

    cube dna(nstate, nstate, ltimes);
    dna.zeros();
    
    for (int t=0; t < ltimes; ++t) {
	for (int i=0; i < nstate; ++i) {
	    if (nrisk.at(t, i) != 0) {
		for (int j=0; j< nstate; ++j) {
		    dna.at(i, j, t) = nev.at(i, j, t) / nrisk.at(t, i);
		}
	    }
	}
	
	mat tmp(dna.slice(t).begin(), nstate, nstate, false);
	vec d = sum(tmp, 1);
    	tmp.diag() = -d;

    }
    
    return dna;

}


/* 
lai and ying modification of the Nelson Aalen estimator.
we have an extra argument: a matrix similar to nrisk, with 0s and 1s
that says where we should compute
*/
cube deltaNA_LY(const cube & nev, const mat & nrisk, const mat & which_compute,
		int nstate, int ltimes) {

    cube dna(nstate, nstate, ltimes);
    dna.zeros();
    
    for (int t=0; t < ltimes; ++t) {
	for (int i=0; i < nstate; ++i) {
	    if (nrisk.at(t, i) != 0) {
		for (int j=0; j< nstate; ++j) {
		    dna.at(i, j, t) = which_compute.at(t, i) * nev.at(i, j, t) / nrisk.at(t, i);
		}
	    }
	}
	
	mat tmp(dna.slice(t).begin(), nstate, nstate, false);
	vec d = sum(tmp, 1);
    	tmp.diag() = -d;

    }
    
    return dna;

}
    

mat cov_dna(const mat & nev, const vec & nrisk, int d, int D) {

    mat the_cov(D, D);
    the_cov.zeros();

    uvec from(D); from.zeros();
    uvec to(D); to.zeros();

    // construct vectors that store the k and l; m, n for each
    // indice of the final covariance matrix
    for (int i = 0; i < d; ++i) {
	for (int j = 0; j < d; ++j) {
	    from[j + i * d] = j;
	    to[j + i * d] = i;
	}
    }
    
    vec sum_nev = sum(nev, 1);
    vec pow_nrisk = pow(nrisk, -3);
    
    for (int j = 0; j < D; ++j) {
	for (int i = 0; i < D; ++i) {

	    if (nrisk[from[i]] != 0) {
		
		int cond = 1 * (from[i] == to[i] && from[i] == from[j] && from[j] == to[j]) +
		    2 * (from[i] == to[i] && from[i] == from[j] && from[i] != to[j]) +
		    4 * (from[i] == from[j] && from[i] != to[i] && from[i] != to[j] && to[i] == to[j]) +
		    8 * (from[i] == from[j] && from[i] != to[i] && from[i] != to[j] && to[i] != to[j]) +
		    // add a condition for the symmetric of cond 2
		    16 * (from[j] == to[j] && from[i] == from[j] && from[i] != to[i]);
		
		switch(cond) {
		case 1:
		    the_cov(i, j) = (nrisk[from[i]] - sum_nev[from[i]]) *
			sum_nev[from[i]] * pow_nrisk[from[i]];
		    break;
		case 2: 
		    the_cov(i, j) = -(nrisk[from[i]] - sum_nev[from[i]]) *
			nev(from[i], to[j]) * pow_nrisk[from[i]];
		    break;
		case 4:
		    the_cov(i, j) = (nrisk[from[i]] - nev(from[i], to[i])) *
			nev(from[i], to[j]) * pow_nrisk[from[i]];
		    break;
		case 8:
		    the_cov(i, j) =  -nev(from[i], to[i]) *
			nev(from[i], to[j]) * pow_nrisk[from[i]];
		    break;
		case 16:
		    the_cov(i, j) = -(nrisk[from[i]] - sum_nev[from[i]]) *
			nev(from[i], to[i]) * pow_nrisk[from[i]];
		    break;
		default:
		    the_cov(i, j) = 0;
		}

	    }
	}
    }

    // Alternative algorithm, just fill the lower triamgular part of the matrix. Then
    // "symmetrize the matrix
    // the_cov = the_cov + the_cov.t();
    // the_cov.diag() /= 2;

//    Rcpp::Rcout << "cov_dna is" << std::endl << the_cov << std::endl;
//    Rcpp::Rcout << "cov_dna[2,5] is" << std::endl << the_cov(1, 4) << std::endl;

    return the_cov;
}
		    
