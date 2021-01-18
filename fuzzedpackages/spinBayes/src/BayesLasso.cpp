#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;

// [[Rcpp::export()]]
Rcpp::List BayesLasso (arma::mat xx, arma::vec y, arma::mat W, arma::mat Ex, unsigned int s, int maxSteps, double hatM, arma::vec hatRStar, arma::vec hatAlpha, arma::vec hatZeta, double invSigM0, arma::vec hatInvTauSqStar, arma::mat invSigAlpha0, arma::vec hatInvTauSqZeta, double hatLambdaSqStar, double hatLambdaSqZeta, double hatSigmaSq, double a0, double b0, double aStar, double bStar, double alpha, double gamma, int progress)
{
	unsigned int n = xx.n_rows, clc = W.n_cols, pR = hatRStar.n_elem, pZ = hatZeta.n_elem;
	arma::mat gsRStar(maxSteps, pR),
			gsAlpha(maxSteps, clc),
			gsZeta(maxSteps, pZ),
			gsInvTauSqZeta(maxSteps, pZ),
			gsInvTauSqStar(maxSteps, pR);
		
	arma::vec gsM(maxSteps),
			gsLambdaZeta(maxSteps),
			gsLambdaStar(maxSteps),
			gsSigmaSq(maxSteps);
		
	arma::mat Br = xx.cols(1, xx.n_cols-1);
	arma::mat tBrBr = Br.t()*Br, tWW = W.t()*W, tExEx = Ex.t()*Ex;
	arma::vec tExExDiag = tExEx.diag(), tBrBrDiag = tBrBr.diag();
	
	arma::mat varAlpha;
	arma::vec res, meanAlpha, tRsRs, muInvTauSqStar, muInvTauSqZeta; // mu_m, mu_alpha, 
	double tempS, tempZ, varM, meanM, meanRs, varRs, meanZ, varZ, ExjtRes, lInvTauSqStar, lInvTauSqZeta;
	
	for (int k = 0; k < maxSteps; k++) {
		if(progress != 0 && k % progress == 0){
			Rcpp::Rcout << "Iteration: " << k << std::endl;
		}
		
		// m|y, r.star
		varM = 1/(n/hatSigmaSq + invSigM0);
		// mu_m = W * hatAlpha + Br * hatRStar;
		res = y - (W * hatAlpha + Ex * hatZeta + Br * hatRStar);
		meanM = varM * arma::accu(res/hatSigmaSq);
		hatM = R::rnorm(meanM, varM);
		res -= hatM;
		gsM(k) = hatM;
		
		// alpha|y, m, r0, r.star
		varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
		// mu_alpha = hatM + Br * hatRStar;
		res += W * hatAlpha;
		meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();
		
		for(unsigned int j=0; j<pZ; j++){
			tempZ = 1/(tExExDiag(j) + hatInvTauSqZeta(j));
			varZ = hatSigmaSq * tempZ;
			res += Ex.col(j) * hatZeta(j);
			ExjtRes = arma::as_scalar(Ex.col(j).t() * res);
			meanZ = tempZ * ExjtRes;
			hatZeta(j) = R::rnorm(meanZ, sqrt(varZ));
			res -= Ex.col(j) * hatZeta(j);
		}
		
		for(unsigned int j=0; j<pR; j++){
			tempS = 1/(tBrBrDiag(j) + hatInvTauSqStar(j));
			varRs = hatSigmaSq * tempS;
			res += Br.col(j) * hatRStar(j);
			meanRs = arma::as_scalar(tempS * Br.col(j).t() * res);
			hatRStar(j) = R::rnorm(meanRs, sqrt(varRs));
			res -= Br.col(j) * hatRStar(j);
		}

		gsZeta.row(k) = hatZeta.t();
		gsRStar.row(k) = hatRStar.t();
		
		// sigma.sq|
		double shapeSig = alpha + (n+pR+pZ)/2;
		// residual = res - Br * hatRStar;
		double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatZeta) % hatInvTauSqZeta) + 
									arma::accu(square(hatRStar) % hatInvTauSqStar));
		hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
		gsSigmaSq(k) = hatSigmaSq;
		
		// invTAUsq.e|
		lInvTauSqZeta = hatLambdaSqZeta;
		muInvTauSqZeta = sqrt(hatLambdaSqZeta * hatSigmaSq / square(hatZeta));
		for(unsigned int j = 0; j < pZ; j++){
			hatInvTauSqZeta(j) = rinvgaussian(muInvTauSqZeta(j), lInvTauSqZeta);
		}
		gsInvTauSqZeta.row(k) = hatInvTauSqZeta.t();
		
		// invTAUsq.star|lambda.star, r.star
		lInvTauSqStar = hatLambdaSqStar;
		tRsRs = square(hatRStar);
		muInvTauSqStar = sqrt(hatLambdaSqStar * hatSigmaSq / tRsRs);
		for(unsigned int j = 0; j<pR; j++){
			hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
		}
		gsInvTauSqStar.row(k) = hatInvTauSqStar.t();
		
		// lambda.zeta|
		double shapeZ = a0 + pZ;
		double rateZ = b0 + arma::accu(1/hatInvTauSqZeta)/2;
		hatLambdaSqZeta = R::rgamma(shapeZ, 1/rateZ);
		gsLambdaZeta(k) = hatLambdaSqZeta;
		
		// lambda.star|invTAUsq.star
		double shapeS = aStar + pR;
		double rateS = bStar + arma::accu(1/hatInvTauSqStar)/2;
		hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
		gsLambdaStar(k) = hatLambdaSqStar;
	}
	
	return Rcpp::List::create(Rcpp::Named("GS.m") = gsM,
							Rcpp::Named("GS.alpha") = gsAlpha,
							Rcpp::Named("GS.zeta") = gsZeta,
							Rcpp::Named("GS.r0") = gsRStar,
							Rcpp::Named("GS.invTAUsq.zeta") = gsInvTauSqZeta,
							Rcpp::Named("GS.invTAUsq.star") = gsInvTauSqStar,
							Rcpp::Named("GS.lambda.sq.zeta") = gsLambdaZeta,
							Rcpp::Named("GS.lambda.sq.star") = gsLambdaStar,
							Rcpp::Named("GS.sigma.sq") = gsSigmaSq);
}
