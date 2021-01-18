#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

Rcpp::List starmaCPP(const arma::mat& data, const Rcpp::List& wlist, arma::mat arMat, arma::mat maMat, int iterate);
Rcpp::List kf(const arma::mat& data, const Rcpp::List& wlist, arma::mat eps, arma::mat arInd, arma::mat maInd);
arma::mat residuals(const arma::mat& data, const Rcpp::List& wlist, const Rcpp::List& model);
double loglik(const arma::mat& data, const Rcpp::List& model);

/* STARMA(p, q, sigma2) parameter estimation function based on a Kalman filter.

   Args:
    - data: A matrix object containing the space time process observations
	- wlist: The list of weight matrices, first one being identity
	- arMat: A 1/0 matrix determining if 'row'-th tlag, 'col'-th slag AR parameter should be estimated
	- maMat: A 1/0 matrix determining if 'row'-th tlag, 'col'-th slag MA parameter should be estimated
	
	Returns:
	A list containing the following elements:
	- $phi: The estimated AR parameters
	- $phi_sd: The estimated standard error of the AR parameters
	- $theta: The estimated MA parameters
	- $theta_sd: The estimated standard error of the MA parameters
	- $sigma2: The estimated white noise variance matrix
   
   The algorithm has been written thanks to the article 'Study on Kalman filter in time series analysis' 
   written by Tomáš Cipra & I. Motyková (1987).
   
   The state space system considered is the following:
		
		ksi[t+1] = ksi[t]
		y[t] = H[t]' * ksi[t] + epsilon[t]
		
   Where ksi is the vector of the parameters,
		 y the vector of observations,
		 H the matrix of past observations,
		 epsilon the vector of innovation.
*/

// [[Rcpp::export]]
Rcpp::List starmaCPP(const arma::mat& data, const Rcpp::List& wlist, arma::mat arMat, arma::mat maMat, int iterate) {

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Create matrices containing the indices of the parameters to estimate
	
	int ar = arma::accu(arMat != 0);	// Nb of AR parameters to estimate
	int ma = arma::accu(maMat != 0);	// Nb of MA parameters to estimate

	// arInd contains the indices of the AR parameters to be estimated
	arma::mat arInd = arma::mat(2, ar);
	int iter = 0;
	for (unsigned int tlag = 0; tlag < arMat.n_rows; tlag++) {
		for (unsigned int slag = 0; slag < arMat.n_cols; slag++) {
			if (arMat(tlag, slag)) {
				arInd(0, iter) = tlag;
				arInd(1, iter) = slag;
				iter++;
			}
		}
	}
	
	// maInd contains the indices of the MA parameters to be estimated
	arma::mat maInd = arma::mat(2, ma);
	iter = 0;
	for (unsigned int tlag = 0; tlag < maMat.n_rows; tlag++) {
		for (unsigned int slag = 0; slag < maMat.n_cols; slag++) {
			if (maMat(tlag, slag)) {
				maInd(0, iter) = tlag;
				maInd(1, iter) = slag;
				iter++;
			}
		}
	}
	
	int p = (arInd.n_elem) ? arma::max(arInd.row(0)) + 1 : 0;	// Max tlag of AR part
	int q = (maInd.n_elem) ? arma::max(maInd.row(0)) + 1 : 0;	// Max tlag of MA part
	int r = (p > q) ? p : q;	// max(p, q)

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Run the Kalman filter.
	// If estimating MA parameters, rerun the Kalman filter, using the previously estimated parameters
	// to get better estimates of the early residuals.
	
	// First iteration
	arma::mat eps = data;	// Estimated residuals
	Rcpp::List model = kf(data, wlist, eps, arInd, maInd);
	
	if (ma) {
		while (iterate--) {
			// Second iteration
			eps.rows(0, r) = residuals(data.rows(0, r), wlist, model);
			model = kf(data, wlist, eps, arInd, maInd);
		}
	}
	
	model["residuals"] = residuals(data, wlist, model);
	model["loglik"] = loglik(data, model);
	model["bic"] = -2 * Rcpp::as<double>(model["loglik"]) + log((double)data.n_elem) * (ar + ma);
	
	return model;

}



Rcpp::List kf(const arma::mat& data, const Rcpp::List& wlist, arma::mat eps, arma::mat arInd, arma::mat maInd) {

	int ar = arInd.n_cols;	// Nb of AR parameters to be estimated
	int ma = maInd.n_cols;	// Nb of MA parameters to be estimated
	int dim = ar + ma;	// Total nb of parameters to be estimated
	int p = (arInd.n_elem) ? arma::max(arInd.row(0)) + 1 : 0;	// Max tlag of AR part
	int l = (arInd.n_elem) ? arma::max(arInd.row(1)) + 1 : 0;	// Max slag of AR part
	int q = (maInd.n_elem) ? arma::max(maInd.row(0)) + 1 : 0;	// Max tlag of MA part
	int m = (maInd.n_elem) ? arma::max(maInd.row(1)) + 1 : 0;	// Max slag of MA part
	int r = (p > q) ? p : q;	// max(p, q)

	// Initialize the Kalman filter
	arma::mat H = arma::mat(dim, data.n_cols);	// Observation matrix
	arma::colvec ksi = arma::zeros(dim);	// Unconditional mean of the prediction set to 0
	arma::mat V = 1e5 * arma::eye(dim, dim);	// V = sigma2^-1 * P, see (Cipra & Motyková, 1987)
	arma::mat sigma2 = 1e-5 * arma::eye(data.n_cols, data.n_cols);	// Variance of the white noise set low
	
	arma::mat weights = arma::mat(data.n_cols, data.n_cols);
	arma::mat Nm1 = arma::mat(data.n_cols, data.n_cols);
	arma::vec nu = arma::vec(data.n_cols);
		
	// Run the filter
	for (unsigned int t = r; t < data.n_rows; t++) {
	
		// Update the observation matrix
		// Fill for 'phi', AR parameters
		for (int it = 0; it < ar; it++) {
			weights = Rcpp::as<arma::mat>(wlist(arInd(1, it)));
			H.row(it) = data.row(t - 1 - arInd(0, it)) * weights.t();
		}
		// Fill for 'theta', MA parameters
		for (int it = ar; it < dim; it++) {
			weights = Rcpp::as<arma::mat>(wlist(maInd(1, it - ar)));
			H.row(it) = eps.row(t - 1 - maInd(0, it - ar)) * weights.t();
		}
	
		Nm1 = arma::inv(H.t() * V * H + arma::eye(data.n_cols, data.n_cols));
		nu = (data.row(t)).t() - H.t() * ksi;
	
		// Prediction & update equations all-in-ones
		ksi += V * H * Nm1 * nu;
		V -= V * H * Nm1 * H.t() * V;
		sigma2 = ((t + 1 - r) * sigma2 + nu * nu.t()) / (t + 2 - r);
		
		// Estimate the residual
		eps.row(t) = data.row(t) - ksi.t() * H;
		
	}
	
	// Get estimated sd of the parameters
	arma::vec sd = arma::sqrt(arma::trace(sigma2) * V.diag() / data.n_cols);
	//arma::vec sd = arma::sqrt(sigma2 * V.diag());
	
	// Rename and reshape 
	arma::mat phi = arma::zeros(p, l);
	arma::mat phi_sd = arma::zeros(p, l);
	arma::mat theta = arma::zeros(q, m);
	arma::mat theta_sd = arma::zeros(q, m);
	
	for (int it = 0; it < ar; it++) {
		phi(arInd(0, it), arInd(1, it)) = ksi[it];
		phi_sd(arInd(0, it), arInd(1, it)) = sd[it];
	}
	
	for (int it = ar; it < dim; it++) {
		theta(maInd(0, it - ar), maInd(1, it - ar)) = ksi[it];
		theta_sd(maInd(0, it - ar), maInd(1, it - ar)) = sd[it];
	}
	
	return Rcpp::List::create(Rcpp::Named("phi")=phi,
							  Rcpp::Named("phi_sd")=phi_sd,
							  Rcpp::Named("theta")=theta,
							  Rcpp::Named("theta_sd")=theta_sd,
							  Rcpp::Named("sigma2")=sigma2);
	
}



/* Iteratively computes the residuals of a given STARMA(p, q) model.

   Args:
    - data: A matrix object containing the space time process observations
	- wlist: The list of weight matrices, first one being identity
	- model: The parameters of the model. This argument should be the output of the 'starmaCPP' function
	
   Returns:
    The residuals of the STARMA(p, q) model.
*/

arma::mat residuals(const arma::mat& data, const Rcpp::List& wlist, const Rcpp::List& model) {
	
	arma::mat phi = Rcpp::as<arma::mat>(model(0));
	arma::mat theta = Rcpp::as<arma::mat>(model(2));
	arma::mat eps = arma::mat(data);	// Will hold the residuals
	arma::mat weights = arma::mat(data.n_cols, data.n_cols);
	int tlim;
	
	// Remove AR part iteratively
	for (unsigned int t = 0; t < data.n_rows; t++) {
		tlim = (t < phi.n_rows) ? t : phi.n_rows;
		for (int tlag = 0; tlag < tlim; tlag++) {
			for (unsigned int slag = 0; slag < phi.n_cols; slag++) {
				weights = Rcpp::as<arma::mat>(wlist(slag));
				eps.row(t) -= phi(tlag, slag) * data.row(t - tlag - 1) * weights.t();
			}
		}
	}
	
	// Remove MA part iteratively
	for (unsigned int t = 0; t < data.n_rows; t++) {
		tlim = (t < theta.n_rows) ? t : theta.n_rows;
		for (int tlag = 0; tlag < tlim; tlag++) {
			for (unsigned int slag = 0; slag < theta.n_cols; slag++) {
				weights = Rcpp::as<arma::mat>(wlist(slag));
				eps.row(t) -= theta(tlag, slag) * eps.row(t - tlag - 1) * weights.t();
			}
		}
	}

	return eps;
	
}



/* Computes the conditional log likelihood of a given STARMA(p, q, sigma2) model.
   
   Args:
    - data: A matrix object containing the space time process observations
	- model: The parameters of the model ; this argument should be the output of the 'starmaCPP' function
	
   Returns:
    The conditional log likelihood of the STARMA(p, q, sigma2) model.
*/

double loglik(const arma::mat& data, const Rcpp::List& model) {

	arma::mat residuals = Rcpp::as<arma::mat>(model(5));
	//arma::mat sigma2 = Rcpp::as<arma::mat>(model(4));		// Cannot use sigma2 variance matrix, as it would skyrocket the nb of parameters to data.n_cols*data._ncols
	double sigma2 = arma::trace(Rcpp::as<arma::mat>(model(4))) / data.n_cols;	
	double L = data.n_cols * data.n_rows * (log(2*arma::datum::pi) + log(sigma2));
	
	for (unsigned int t = 0; t < data.n_rows; t++) {
		//L += arma::as_scalar(residuals.row(t) * inv(sigma2) * (residuals.row(t)).t());
		L += arma::as_scalar(residuals.row(t) * (1/sigma2) * (residuals.row(t)).t());
	}
	
	return -L/2;
	
}
