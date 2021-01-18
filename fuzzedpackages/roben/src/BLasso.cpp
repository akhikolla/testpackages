#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BLasso (arma::mat xx, arma::vec y, arma::mat W, int maxSteps, arma::vec hatBeta, arma::vec hatAlpha, arma::vec hatInvTauSq, arma::mat invSigAlpha0, double hatLambdaSqStar, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, int progress)
{
	unsigned int n = xx.n_rows, s = xx.n_cols, clc = W.n_cols;
	arma::mat gsAlpha(maxSteps, clc),
			gsBeta(maxSteps, s),
			gsInvTauSq(maxSteps, s);
		
	arma::vec gsLambdaStar(maxSteps),
			gsSigmaSq(maxSteps),
			gsMSE(maxSteps);

	arma::mat tWW = W.t()*W, varAlpha, matBeta;
	arma::vec res, meanAlpha, tRsRs, muInvTauSq;
	double tempS, varRs, meanRs, lInvTauSq;
	
	arma::vec tBrBrDiag = sum(square(xx), 0).t();
	
	for (int k = 0; k < maxSteps; k++) {
		// alpha|
		varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
		res = y - xx * hatBeta;
		meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();

		for(unsigned int j=0; j<s; j++){		
			tempS = 1/(tBrBrDiag(j) + hatInvTauSq(j));
			varRs = hatSigmaSq * tempS;
			res += xx.col(j) * hatBeta(j);
			meanRs = arma::as_scalar(tempS * xx.col(j).t() * res);
			hatBeta(j) = R::rnorm(meanRs, std::sqrt(varRs));
			res -= xx.col(j) * hatBeta(j);
		}
		gsBeta.row(k) = hatBeta.t();
		
		// sigma.sq|
		double shapeSig = alpha + (n+s)/2;
		double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatBeta) % hatInvTauSq));
		hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
		gsSigmaSq(k) = hatSigmaSq;
		
		
		// invTAUsq.star|
		lInvTauSq = hatLambdaSqStar;
		tRsRs = arma::square(hatBeta);
		muInvTauSq = arma::sqrt(hatLambdaSqStar * hatSigmaSq / tRsRs);
		for(unsigned int j = 0; j<s; j++){
			hatInvTauSq(j) = rinvgaussian(muInvTauSq(j), lInvTauSq);
		}
		gsInvTauSq.row(k) = hatInvTauSq.t();
		
		
		// lambda.star|
		double shapeS = aStar + s;
		double rateS = bStar + arma::accu(1/hatInvTauSq)/2;
		hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
		gsLambdaStar(k) = hatLambdaSqStar;
		
		gsMSE(k) = arma::mean(arma::square(res));
		if(k % 100 == 0){
			Rcpp::checkUserInterrupt();
		}
		if(progress != 0 && k % progress == 0){
			Rcpp::Rcout << "\nIter." << k << "  mse: " << gsMSE(k) << std::endl;
		}
	}
	
	return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
							Rcpp::Named("GS.beta") = gsBeta,
							Rcpp::Named("GS.invTAUsq") = gsInvTauSq,
							Rcpp::Named("GS.lambda.sq") = gsLambdaStar,
							Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
							Rcpp::Named("GS.mse") = gsMSE);
}
