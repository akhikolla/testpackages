#include "generate_norm.h"

using namespace Rcpp;

SEXP rmvnorm(SEXP nsims, SEXP mus, SEXP Vs, SEXP methods){

	int nsim = as<int>(nsims);
	
	NumericVector mur(mus);
	arma::vec mu(mur.begin(), mur.size(), false);
	
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
	
	RNGScope scope;
	NumericVector Zr = NumericVector(rnorm(Vr.nrow() * nsim, 0, 1));
	arma::mat Z(Zr.begin(), Vr.nrow(), nsim);
	
	return(wrap(repmat(mu, 1, nsim) + decomp_V * Z));
}

SEXP condnorm_par(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP coeffs, SEXP Xs, SEXP Xps, SEXP methods){

	NumericVector yr(ys);
	arma::colvec y(yr.begin(), yr.size(), false);

	NumericMatrix Vr(Vs);
	arma::mat V(Vr.begin(), Vr.nrow(), Vr.ncol(), false);

	NumericMatrix Vpr(Vps);
	arma::mat Vp(Vpr.begin(), Vpr.nrow(), Vpr.ncol(), false);

	NumericMatrix Vopr(Vops);
	arma::mat Vop(Vopr.begin(), Vopr.nrow(), Vopr.ncol(), false);

	NumericVector coeffr(coeffs);
	arma::colvec coeff(coeffr.begin(), coeffr.size(), false);

	NumericMatrix Xr(Xs);
	arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);

	NumericMatrix Xpr(Xps);
	arma::mat Xp(Xpr.begin(), Xpr.nrow(), Xpr.ncol(), false);

	int method = as<int>(methods);

	arma::mat ViVop = solve(V, Vop);

	arma::mat mc = Xp * coeff + trans(ViVop) * (y - X * coeff);
	arma::mat Vc = Vp - trans(Vop) * ViVop;

	arma::mat decomp_Vc = arma::mat(Vp.n_rows, Vp.n_rows);

	if(method == 1){
	arma::vec eigval = arma::vec(Vp.n_rows);
	arma::mat eigvec = arma::mat(Vp.n_rows, Vp.n_rows);

	//compute eigen values and vectors of Vc
	eig_sym(eigval, eigvec, Vc);

	for(int i = 0; i < eigval.n_rows; i++)
	{
	if(eigval(i) < 0)
	{
		eigval(i) = 0;		
	}
	}

	decomp_Vc = eigvec * diagmat(sqrt(eigval));
	}
	else if(method == 2){
	decomp_Vc = trans(chol(Vc));
	}
	else{
	arma::mat U = arma::mat(Vp.n_rows, Vp.n_rows);
	arma::mat U2 = arma::mat(Vp.n_rows, Vp.n_rows);
	arma::vec sv = arma::vec(Vp.n_rows);

	svd(U, sv, U2, Vc);

	decomp_Vc = U * diagmat(sqrt(sv)) * trans(U2);
	}

	/*//RNGScope scope;
	NumericVector Zr = NumericVector(rnorm(Vr.nrow() * nsim, 0, 1));
	arma::mat Z(Zr.begin(), Vr.nrow(), nsim);
	*/

	return Rcpp::List::create(Rcpp::Named("mc") = mc,
					  Rcpp::Named("decomp.Vc")= decomp_Vc
					  );
}
