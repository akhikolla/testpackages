#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BGL (arma::mat xx, arma::vec y, arma::mat W, unsigned int s, unsigned int q, int maxSteps, arma::vec hatBeta, arma::vec hatAlpha, arma::vec hatInvTauSq, arma::mat invSigAlpha0, double hatLambdaSqStar, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, int progress)
{
	unsigned int L = q, n = xx.n_rows, clc = W.n_cols;
	arma::mat gsAlpha(maxSteps, clc),
			gsBeta(maxSteps, s*q),
			gsInvTauSq(maxSteps, s);
		
	arma::vec gsLambdaStar(maxSteps),
			gsSigmaSq(maxSteps),
			gsMSE(maxSteps);

	arma::mat Xr, tWW = W.t()*W, varAlpha, varRs, tempS, matBeta;
	arma::vec res, meanAlpha, meanRs, BrjtRes, tRsRs, repInvTau, muInvTauSq;
	double lInvTauSq;
	
	std::vector<arma::mat> tBrBr(s);
	for(unsigned int j=0; j<s; j++){
		Xr = xx.cols((j*L), (j*L+L-1));
		tBrBr[j] = Xr.t()*Xr;
	}
	
	for (int k = 0; k < maxSteps; k++) {
		// alpha|
		varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
		res = y - xx * hatBeta;
		meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();
		
		// beta|
		for(unsigned int j=0; j<s; j++){		
			tempS = tBrBr[j];
			tempS.diag() += hatInvTauSq(j);
			tempS = arma::inv(tempS);
			varRs = hatSigmaSq * tempS;
			res += xx.cols((j*L), (j*L+L-1)) * hatBeta.subvec((j*L), (j*L+L-1));
			BrjtRes = xx.cols((j*L), (j*L+L-1)).t() * res;
			meanRs = tempS * BrjtRes;
			hatBeta.subvec((j*L), (j*L+L-1)) = mvrnormCpp(meanRs, varRs);
			res -= xx.cols((j*L), (j*L+L-1)) * hatBeta.subvec((j*L), (j*L+L-1));
		}
		gsBeta.row(k) = hatBeta.t();
		
		// sigma.sq|
		double shapeSig = alpha + (n+s*L)/2;
		repInvTau = arma::vectorise(arma::repelem(hatInvTauSq.t(), L, 1), 0);
		double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatBeta) % repInvTau));
		hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
		gsSigmaSq(k) = hatSigmaSq;
		
		
		// invTAUsq.star|
		lInvTauSq = L * hatLambdaSqStar;
		matBeta = arma::reshape(hatBeta, L, s);
		tRsRs = sum(square(matBeta), 0).t();
		muInvTauSq = sqrt(L * hatLambdaSqStar * hatSigmaSq / tRsRs);
		for(unsigned int j = 0; j<s; j++){
			hatInvTauSq(j) = rinvgaussian(muInvTauSq(j), lInvTauSq);
		}
		gsInvTauSq.row(k) = hatInvTauSq.t();
		
		
		// lambda.star|
		double shapeS = aStar + s*(L+1)/2;
		double rateS = bStar + L*arma::accu(1/hatInvTauSq)/2;
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
