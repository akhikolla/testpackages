#include "GaussianProcess.h"

GaussianProcess::GaussianProcess(int Inputs, int Outputs, mat& Xdata,
		vec& ydata, CovarianceFunction& cf) :
    ForwardModel(Inputs, Outputs), covFunc(cf), Locations(Xdata), 
               Observations(ydata)  {
	

}

GaussianProcess::~GaussianProcess() {
}

//void GaussianProcess::train(mat& C, const mat& X) const
//{
//	mat A;
//	A.set_size(C.n_rows, C.n_cols, false);
//	covFunc->computeSymmetric(A, X);
//	inv(A, C);
//}

void GaussianProcess::makePredictions(vec& Mean, vec& Variance,
		const mat& Xpred) const {
	makePredictions(Mean, Variance, Xpred, covFunc);
}

void GaussianProcess::makePredictions(vec& Mean, vec& Variance,
		const mat& Xpred, CovarianceFunction &cf) const {
	
	

	mat Sigma(Observations.size(), Observations.size());
	mat Cpred(Locations.n_rows, Xpred.n_rows);

	cf.computeCovariance(Cpred, Locations, Xpred); // k = k(X,x*)
	covFunc.computeSymmetric(Sigma, Locations); // K = K(X,X)

	// vec alpha = ls_solve(Sigma, Observations);                       // a = K^{-1} * y
	vec alpha = arma::solve(inv(Sigma), Observations);
	Mean = Cpred.t() * alpha; // mu* = k' * K^{-1} * y

	vec variancePred(Xpred.n_rows);
	cf.computeDiagonal(variancePred, Xpred); // k* = K(x*,x*)
	// mat v = ls_solve(computeCholesky(Sigma).transpose(), Cpred);     // v = K^{-1} * k
	mat v = solve(inv(computeCholesky(Sigma).t()), Cpred);
	Variance = variancePred - sum(v % v); // diag( k* - k'*K^{-1}*k )
}

void GaussianProcess::makePredictions(vec& Mean, vec& Variance,
		const mat& Xpred, const mat& C) const {

	
	

	mat Sigma, cholSigma;
	cholSigma.set_size(Observations.size(), Observations.size());

	Sigma.set_size(Observations.size(), Observations.size());
//	cholSigma.set_size(Observations.n_rows, Observations.n_rows);

	covFunc.computeSymmetric(Sigma, Locations);

	mat Cpred(Locations.n_rows, Xpred.n_rows);

	covFunc.computeCovariance(Cpred, Locations, Xpred);

	//mat alpha = C * Observations;

	covFunc.computeSymmetric(Sigma, Locations);

	cholSigma = computeCholesky(Sigma);

	// vec alpha = ls_solve_chol(cholSigma, Observations);
	vec alpha = solve(inv(Sigma), Observations);

	Mean = Cpred.t() * alpha;

	vec sigsq(Xpred.n_rows);

	covFunc.computeDiagonal(sigsq, Xpred);

	Variance = sum((Cpred.t() * C) % Cpred.t(), 0);
	Variance = sigsq + Variance;

}

double GaussianProcess::loglikelihood() const {

	mat Sigma(Observations.size(), Observations.size());
	mat cholSigma(Observations.size(), Observations.size());

	covFunc.computeSymmetric(Sigma, Locations);

	cholSigma = computeCholesky(Sigma);

	mat invSigma = computeInverseFromCholesky(Sigma);
	vec alpha = invSigma * Observations;

//	vec alpha = ls_solve(Sigma, Observations);

	double out1 = 0.5 * dot(Observations, alpha);

	double out2 = arma::accu(arma::log(arma::diagvec(cholSigma)));

	return out1 + out2 + 0.5 * Observations.n_elem * log(2 * arma::datum::pi);

}

vec GaussianProcess::getParametersVector() const {
	return covFunc.getParameters();
}

void GaussianProcess::setParametersVector(const vec p) {
	covFunc.setParameters(p);
}

double GaussianProcess::objective() const {
	return loglikelihood();
}

/*
 vec GaussianProcess::gradient() const
 {
 vec grads(covFunc.getNumberParameters());


 mat Sigma(Observations.size(), Observations.size());
 mat cholSigma(Observations.size(), Observations.size());

 covFunc.computeSymmetric(Sigma, Locations);

 cholSigma = computeCholesky(Sigma);

 vec alpha = ls_solve(Sigma, Observations);
 //	mat W = ls_solve(cholSigma, ls_solve(cholSigma.transpose(), eye(Observations.size()))) - outer_product(alpha, alpha, false);
 mat W = (inv(Sigma) - outer_product(alpha, alpha, false));

 mat partialDeriv(Observations.size(), Observations.size());

 for(int i = 0; i < covFunc.getNumberParameters(); i++)
 {
 covFunc.getParameterPartialDerivative(partialDeriv, i, Locations);
 //grads(i) = sum(sum(elem_mult(W, partialDeriv))) / 2;
 grads(i) = elem_mult_sum(W, partialDeriv) / 2;
 // official - but slower		grads(i) = sum(diag(W * partialDeriv)) / 2;
 }
 return grads;
 }
 */

vec GaussianProcess::gradient() const {
	vec grads(covFunc.getNumberParameters());

	mat Sigma(Observations.n_elem, Observations.n_elem);
	mat cholSigma(Observations.n_elem, Observations.n_elem);

	covFunc.computeSymmetric(Sigma, Locations);
	cholSigma = computeCholesky(Sigma);
	mat invSigma = computeInverseFromCholesky(Sigma);
	vec alpha = invSigma * Observations;

//	vec alpha = ls_solve(Sigma, Observations);
	mat W = (invSigma - alpha * alpha.t());

	mat partialDeriv(Observations.size(), Observations.size());

	for (unsigned int i = 0; i < covFunc.getNumberParameters(); i++) {
		covFunc.getParameterPartialDerivative(partialDeriv, i, Locations);
		//grads(i) = sum(sum(elem_mult(W, partialDeriv))) / 2;
		grads(i) = arma::accu(W % partialDeriv) / 2.0;
// official - but slower		grads(i) = sum(diag(W * partialDeriv)) / 2;
	}
	return grads;
}

vec GaussianProcess::getGradientVector() const {
	return gradient();
}

mat GaussianProcess::computeCholesky(const mat& iM) const {
	mat M = iM; // oops, was i inadvertantly writing to this?
	

	const double ampl = 1.0e-10;
	const int maxAttempts = 10;

	mat cholFactor(M.n_rows, M.n_cols);

	int l = 0;
	bool success = chol(M, cholFactor);
	if (success) {
		return cholFactor;
	} else {
		double noiseFactor = abs(ampl * (trace(M) / double(M.n_rows)));
		while (!success) {
			M = M + (noiseFactor * arma::eye(M.n_rows, M.n_rows));

			if (l > maxAttempts) {
				Rprintf("Unable to compute cholesky decomposition");
				break;
			}
			l++;
			noiseFactor = noiseFactor * 10;
			success = chol(M, cholFactor);
		}
		Rprintf(
				"Matrix not positive definite.  After %d attempts, %f added to the diagonal",
				l, noiseFactor);
	}
	return cholFactor;
}

mat GaussianProcess::computeInverseFromCholesky(const mat& C) const {
	mat cholFactor = computeCholesky(C);
	mat invChol = solve(cholFactor,
			arma::eye(cholFactor.n_rows, cholFactor.n_rows));
	return invChol * invChol.t();
}

