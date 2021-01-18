#include "cov_tools.h"

using namespace Rcpp;

SEXP decomp_cov(SEXP Vs, SEXP methods){

	NumericMatrix Vr(Vs);
	arma::mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);
	
	int method = as<int>(methods);
	
	arma::mat decomp_V = arma::mat(Vr.nrow(), Vr.nrow());
	
	if(method == 1){
		arma::vec eigval = arma::vec(Vr.nrow());
		arma::mat eigvec = arma::mat(Vr.nrow(), Vr.nrow());
		
		//compute eigen values and vectors of V
		eig_sym(eigval, eigvec, V);
		
		for(int i = 0; i < eigval.n_rows; i++)
		{
			if(eigval(i) < 0)
			{
				eigval(i) = 0;		
			}
		}
		
		decomp_V = eigvec * diagmat(sqrt(eigval));
	}
	else if(method == 2){
		decomp_V = trans(chol(V));
	}
	else{
		arma::mat U = arma::mat(Vr.nrow(), Vr.nrow());
		arma::mat U2 = arma::mat(Vr.nrow(), Vr.nrow());
		arma::vec sv = arma::vec(Vr.nrow());
		
		svd(U, sv, U2, V);
		
		decomp_V = U * diagmat(sqrt(sv)) * trans(U2);
    }
	
	return(wrap(decomp_V));
}
