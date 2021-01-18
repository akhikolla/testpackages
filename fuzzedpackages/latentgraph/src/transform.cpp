#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export]]
RcppExport SEXP transformEachObs(SEXP Rj, SEXP Rp, SEXP RR, SEXP RX){
int j = as<int>(Rj);
int p = as<int>(Rp);
int R = as<int>(RR);
arma::mat X = Rcpp::as<arma::mat>(RX);
arma::vec y = X.col(j-1);
arma::mat Z = X;
Z.shed_col(j-1);

arma::vec ytilde(Rf_choose(R,2));
arma::mat Xtilde(Rf_choose(R,2),(p-1));

int i = 0;
for(int k=0; k<(R-1); k++){
	for(int l=(k+1); l<R; l++){
		ytilde[i] = y[k]-y[l];
		Xtilde.row(i) = Z.row(k)-Z.row(l);
		i++;
	}
}

arma::vec finaly(Rf_choose(R,2));
arma::mat finalX(Rf_choose(R,2),(p-1));
arma::vec temp = arma::ones(p-1);
finaly = sign(ytilde);
finalX = Xtilde%abs(ytilde*temp.t());

return Rcpp::List::create(Rcpp::Named("ytilde")=finaly,Rcpp::Named("Xtilde")=finalX);
}
