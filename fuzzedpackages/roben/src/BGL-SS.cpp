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
Rcpp::List BGLPointMass (arma::mat xx, arma::vec y, arma::mat W, unsigned int s, unsigned int q, int maxSteps, arma::vec hatAlpha, arma::vec hatBeta, arma::vec hatInvTauSqStar, arma::mat invSigAlpha0, double hatPiStar, double hatLambdaSqStar, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, double sh1, double sh0, int progress)
{
	unsigned int L = q, n = xx.n_rows, clc = W.n_cols;
	arma::mat gsAlpha(maxSteps, clc),
			gsBeta(maxSteps, s*q),
			gsRstRs(maxSteps, s),
			gsInvTauSqStar(maxSteps, s),
			gsLS(maxSteps, s);
		
	arma::vec gsLambdaStar(maxSteps),
			gsSigmaSq(maxSteps),
			gsPiStar(maxSteps),
			gsMSE(maxSteps);

	arma::mat Br = xx;
	arma::mat tWW = W.t()*W;
	
	arma::mat Xr, varM, varAlpha, varRs, tempS, matRStar;
	arma::vec res, BrjtRes, meanAlpha, meanRs, tRsRs, repInvTauSq, muInvTauSqStar; 
	double lS, t, lInvTauSqStar;
	
	std::vector<arma::mat> tBrBr(s);
	for(unsigned int j=0; j<s; j++){
		Xr = Br.cols((j*L), (j*L+L-1));
		tBrBr[j] = Xr.t()*Xr;
	}

	for (int k = 0; k < maxSteps; k++) {
		// alpha|
		varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
		res = y - Br * hatBeta;
		meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();
		
		// beta|
		for(unsigned int j=0; j<s; j++){		
			tempS = tBrBr[j];
			tempS.diag() += hatInvTauSqStar(j);
			tempS = arma::inv(tempS);
			varRs = hatSigmaSq * tempS;
			res += Br.cols((j*L), (j*L+L-1)) * hatBeta.subvec((j*L), (j*L+L-1));
			BrjtRes = Br.cols((j*L), (j*L+L-1)).t() * res;
			meanRs = tempS * BrjtRes;
			lS = arma::as_scalar(hatPiStar/(hatPiStar+(1-hatPiStar)*std::pow(hatInvTauSqStar(j),(L/2))*std::sqrt(arma::det(tempS))*arma::exp(0.5/hatSigmaSq*(BrjtRes.t()*tempS*BrjtRes))));
			gsLS(k, j) = lS;
			t = R::runif(0, 1);
			if(t<lS){
				hatBeta.subvec((j*L), (j*L+L-1)).zeros();
			}else{
				hatBeta.subvec((j*L), (j*L+L-1)) = mvrnormCpp(meanRs, varRs);
			}
			res -= Br.cols((j*L), (j*L+L-1)) * hatBeta.subvec((j*L), (j*L+L-1));
		}
		gsBeta.row(k) = hatBeta.t();
				
		// sigma.sq|
		double shapeSig = alpha + n/2 + L*arma::accu(tRsRs != 0)/2;
		repInvTauSq = arma::vectorise(arma::repelem(hatInvTauSqStar.t(), L, 1), 0);
		double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatBeta) % repInvTauSq));
		hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
		gsSigmaSq(k) = hatSigmaSq;
		
		// invTAUsq.star|
		lInvTauSqStar = L * hatLambdaSqStar;
		matRStar = arma::reshape(hatBeta, L, s);
		tRsRs = arma::sum(arma::square(matRStar), 0).t();
		muInvTauSqStar = arma::sqrt(L * hatLambdaSqStar * hatSigmaSq / tRsRs);		
		for(unsigned int j = 0; j<s; j++){
			if(tRsRs(j) == 0){
				hatInvTauSqStar(j) = 1/R::rgamma((L+1)/2, 2/lInvTauSqStar);
			}else{
				hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
			}
		}
		gsInvTauSqStar.row(k) = hatInvTauSqStar.t();
		gsRstRs.row(k) = tRsRs.t();
	
		// lambda.star|
		double shapeS = aStar + s*(L+1)/2;
		double rateS = bStar + L*arma::accu(1/hatInvTauSqStar)/2;
		hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
		gsLambdaStar(k) = hatLambdaSqStar;

		
		// pi.star|
		double shape1_s = sh0 + arma::accu(tRsRs == 0);
		double shape2_s = sh1 + arma::accu(tRsRs != 0);
		hatPiStar = R::rbeta(shape1_s, shape2_s);
		gsPiStar(k) = hatPiStar;
		
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
								Rcpp::Named("GS.tRsRs") = gsRstRs,
								Rcpp::Named("GS.invTAUsq") = gsInvTauSqStar,
								Rcpp::Named("GS.pi") = gsPiStar,
								Rcpp::Named("GS.lambda.sq") = gsLambdaStar,
								Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
								Rcpp::Named("GS.lS") = gsLS);
}
