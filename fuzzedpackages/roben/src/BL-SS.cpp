#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include<vector>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BL_SS (arma::mat xx, arma::vec y, arma::mat W, int maxSteps, arma::vec hatAlpha, arma::vec hatBeta, arma::vec hatInvTauSq, arma::mat invSigAlpha0, double hatPi, double hatLambdaSq, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, double sh1, double sh0, int progress)
{
	unsigned int n = xx.n_rows, s = xx.n_cols, clc = W.n_cols;
	arma::mat gsAlpha(maxSteps, clc),
			gsBeta(maxSteps, s),
			gsInvTauSqStar(maxSteps, s),
			gsLS(maxSteps, s);
		
	arma::vec gsLambda(maxSteps),
			gsSigmaSq(maxSteps),
			gsPiStar(maxSteps),
			gsMSE(maxSteps);

	arma::mat tWW = W.t()*W;
	arma::mat varM, varAlpha, matRStar;
	arma::vec tXX = arma::sum(arma::square(xx), 0).t(), res, meanAlpha, muInvTauSqStar; 
	double tempS, meanRs, XjTRes, varRs, lS, t, lInvTauSqStar;

	for (int k = 0; k < maxSteps; k++) {
		// Rcpp::Rcout << "alpha" << std::endl;
		varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
		res = y - xx * hatBeta;
		meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();
		
		// Rcpp::Rcout << "beta" << std::endl;
		for(unsigned int j=0; j<s; j++){
			res += xx.col(j) * hatBeta(j);
			tempS = 1/(tXX(j) + hatInvTauSq(j));
			varRs = hatSigmaSq * tempS;
			XjTRes = arma::as_scalar(xx.col(j).t() * res);
			meanRs = tempS * XjTRes;
			double lS_temp = std::exp(-0.5/hatSigmaSq*tempS*XjTRes*XjTRes)/std::sqrt(hatInvTauSq(j))/std::sqrt(tempS);
			lS = hatPi/(hatPi+(1-hatPi)*lS_temp);
			gsLS(k, j) = lS;
			t = R::runif(0, 1);
			if(t<lS){
				hatBeta(j) = R::rnorm(meanRs, std::sqrt(varRs));
			}else{
				hatBeta(j) = 0;
			}
			res -= xx.col(j) * hatBeta(j);
		}
		gsBeta.row(k) = hatBeta.t();
	
		// Rcpp::Rcout << "invTAUsq.star" << std::endl;
		lInvTauSqStar = hatLambdaSq;
		muInvTauSqStar = std::sqrt(hatLambdaSq * hatSigmaSq) / arma::abs(hatBeta);		
		for(unsigned int j = 0; j<s; j++){
			if(hatBeta(j) == 0){
				hatInvTauSq(j) = 1/R::rexp(2/lInvTauSqStar);
			}else{
				hatInvTauSq(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
			}
		}
		gsInvTauSqStar.row(k) = hatInvTauSq.t();
	
		// lambda.star|
		double rate = bStar + arma::accu(1/hatInvTauSq)/2;
		hatLambdaSq = R::rgamma(aStar+s, 1/rate);
		gsLambda(k) = hatLambdaSq;
		
		// pi.star|
		double shape1 = sh1 + arma::accu(hatBeta != 0);
		double shape2 = sh0 + arma::accu(hatBeta == 0);
		hatPi = R::rbeta(shape1, shape2);
		gsPiStar(k) = hatPi;
		
		// sigma.sq|
		double shapeSig = alpha + n/2 + arma::accu(hatBeta != 0)/2;		
		double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatBeta) % hatInvTauSq));
		hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
		gsSigmaSq(k) = hatSigmaSq;
		
		
		gsMSE(k) = arma::mean(arma::square(res));
		if(k % 100 == 0){
			Rcpp::checkUserInterrupt();
		}
		if(progress != 0 && k % progress == 0){
			Rcpp::Rcout << "\nIter." << k << "  mse: " << gsMSE(k) << std::endl;
		}
	}
	return Rcpp::List::create(	Rcpp::Named("GS.alpha") = gsAlpha,
								Rcpp::Named("GS.beta") = gsBeta,
								Rcpp::Named("GS.invTAUsq") = gsInvTauSqStar,
								Rcpp::Named("GS.pi") = gsPiStar,
								Rcpp::Named("GS.lambda.sq") = gsLambda,
								Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
								Rcpp::Named("GS.lS") = gsLS);
}
