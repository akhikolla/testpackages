#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double stcovCPP(const arma::mat& data, const Rcpp::List& wlist, int slag1, int slag2, int tlag) {

	double out = 0;
	int indlim = data.n_rows - tlag;
	
	arma::mat w1 = Rcpp::as<arma::mat>(wlist(slag1));
	arma::mat w2 = Rcpp::as<arma::mat>(wlist(slag2));
	arma::mat pW = w2.t() * w1;
	
	for (int t = 0; t < indlim; t++) {
		out += trace( pW * trans(data.row(t)) * data.row(t+tlag) );
	}
	
	out /= indlim * data.n_cols;
	
	return out;
}

// [[Rcpp::export]]
arma::mat stacfCPP(const arma::mat& data, const Rcpp::List& wlist, int tlagMax) {

	arma::mat out = arma::mat(tlagMax, wlist.size());
	
	double cov000 = stcovCPP(data, wlist, 0, 0, 0);
	for (int slag = 0; slag < wlist.size(); slag++) {
		double covss0 = stcovCPP(data, wlist, slag, slag, 0);
		for (int tlag = 1; tlag <= tlagMax; tlag++) {
			double covs0t = stcovCPP(data, wlist, slag, 0, tlag);
			out(tlag - 1, slag) = covs0t / sqrt( covss0 * cov000 );
		}
	}
	
	return out;
}

arma::mat stmatCPP_(const arma::mat& data, const Rcpp::List& wlist, int tlag) {
	
	int slagMax = wlist.size();
	arma::mat out = arma::mat(slagMax, slagMax);
	
	for (int i = 0; i < slagMax; i++) {
		for (int j = 0; j < slagMax; j++) {
			out(i, j) = stcovCPP(data, wlist, i, j, tlag);
		}
	}
	
	return out;
	
}

arma::mat stmatCPP(const arma::mat& data, const Rcpp::List& wlist, int tlagMax) {
	
	int slagMax = wlist.size();
	arma::mat slideye = arma::eye(tlagMax, 2*tlagMax - 1);
	arma::mat out = arma::zeros(slagMax * tlagMax, slagMax * tlagMax);

	for (int tlag = 1; tlag < tlagMax; tlag++) {
		out += arma::kron(slideye.submat(0, tlag, tlagMax - 1, tlag + tlagMax - 1), stmatCPP_(data, wlist, tlag));
	}
	
	out += out.t();
	out += arma::kron(arma::eye(tlagMax, tlagMax), stmatCPP_(data, wlist, 0));
	
	return out;
	
}

arma::colvec stvecCPP(const arma::mat& data, const Rcpp::List& wlist, int tlagMax) {
	
	int slagMax = wlist.size();
	arma::colvec out = arma::colvec(slagMax * tlagMax);
	
	for (int tlag = 1; tlag <= tlagMax; tlag++) {
		for (int slag = 0; slag < slagMax; slag++) {
			out[(tlag - 1) * slagMax + slag] = stcovCPP(data, wlist, slag, 0, tlag);
		}
	}
	
	return out;
	
}

// [[Rcpp::export]]
arma::mat stpacfCPP(const arma::mat& data, const Rcpp::List& wlist, int tlagMax) {
	
	int slagMax = wlist.size();
	arma::mat YWmat = stmatCPP(data, wlist, tlagMax);
	arma::colvec YWvec = stvecCPP(data, wlist, tlagMax);
	arma::mat out = arma::mat(tlagMax, slagMax);
	
	for (int tlag = 1; tlag <= tlagMax; tlag++) {
		for (int slag = 0; slag < slagMax; slag++) {
			int index = (tlag - 1) * slagMax + slag;
			arma::colvec sol = solve(YWmat.submat(0, 0, index, index), YWvec.subvec(0, index));
			out(tlag - 1, slag) = sol[index];
		}
	}
	
	return out;
	
}
