#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BSGL (arma::mat xx, arma::vec y, arma::mat W, unsigned int s, unsigned int L, int maxSteps, arma::vec hatAlpha, arma::mat hatBeta, arma::vec hatInvTauSq, arma::mat hatInvGammaSq, arma::mat invSigAlpha0, double hatLambdaSqT, double hatLambdaSqG, double hatSigmaSq, double s1, double s2, double r1, double r2, double a, double b, int progress)
{
	unsigned int n = xx.n_rows, clc = W.n_cols;
	arma::mat gsAlpha(maxSteps, clc),
			gsBeta(maxSteps, s*L),
			gsInvGammaSq(maxSteps, s*L),
			gsInvTauSq(maxSteps, s);
		
	arma::vec gsLambdaSqT(maxSteps),
			gsLambdaSqG(maxSteps),
			gsSigmaSq(maxSteps),
			gsMSE(maxSteps);

	arma::mat tWW = W.t()*W;
	arma::mat Xr, varAlpha, varRs, tempS, repInvTau;
	arma::vec res, meanAlpha, meanRs, muG, tRsRs, muTau;
	
	std::vector<arma::mat> tBrBr(s);
	for(unsigned int j=0; j<s; j++){
		Xr = xx.cols((j*L), (j*L+L-1));
		tBrBr[j] = Xr.t()*Xr;
	}
	
	for (int k = 0; k < maxSteps; k++) {
		// Rcpp::Rcout << "alpha" << std::endl;
		res = y - xx * arma::vectorise(hatBeta);
		varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
		meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();
		
		// Rcpp::Rcout << "beta" << std::endl;
		for(unsigned int j=0; j<s; j++){
			res += xx.cols((j*L), (j*L+L-1)) * hatBeta.col(j);
			
			tempS = tBrBr[j];
			tempS.diag() += (hatInvGammaSq.col(j) + hatInvTauSq(j));
			tempS = arma::inv(tempS);
			
			varRs = hatSigmaSq * tempS;
			
			meanRs = tempS * xx.cols((j*L), (j*L+L-1)).t() * res;
			hatBeta.col(j) = mvrnormCpp(meanRs, varRs);
			res -= xx.cols((j*L), (j*L+L-1)) * hatBeta.col(j);
		}
		gsBeta.row(k) = arma::vectorise(hatBeta).t();
		
		// Rcpp::Rcout << "invGammSq" << std::endl;
		for(unsigned int j = 0; j<s; j++){
			muG = std::sqrt(hatLambdaSqG * hatSigmaSq) / arma::abs(hatBeta.col(j));
			for(unsigned int l = 0; l<L; l++){
				hatInvGammaSq(l,j) = rinvgaussian(muG(l), hatLambdaSqG);
			}
		}
		gsInvGammaSq.row(k) = arma::vectorise(hatInvGammaSq).t();
		
		// Rcpp::Rcout << "invTauSq" << std::endl;
		tRsRs = arma::sum(square(hatBeta), 0).t();
		muTau = arma::sqrt(hatLambdaSqT * hatSigmaSq / tRsRs);
		for(unsigned int j = 0; j<s; j++){
			hatInvTauSq(j) = rinvgaussian(muTau(j), hatLambdaSqT);
		}
		gsInvTauSq.row(k) = hatInvTauSq.t();
		
		// Rcpp::Rcout << "Lambda.Gamma" << std::endl;
		double shapeG = s*L + s1;
		double rateG = r1 + arma::accu(1/hatInvGammaSq)/2;
		hatLambdaSqG = R::rgamma(shapeG, 1/rateG);
		gsLambdaSqG(k) = hatLambdaSqG;
		
		// Rcpp::Rcout << "Lambda.Tau" << std::endl;
		double shapeT = s/2 + s2;
		double rateT = r2 + arma::accu(1/hatInvTauSq)/2;
		hatLambdaSqT = R::rgamma(shapeT, 1/rateT);
		gsLambdaSqT(k) = hatLambdaSqT;
		
		// Rcpp::Rcout << "sigma.sq" << std::endl;
		double shapeSig = a + (n+s*L)/2;
		repInvTau = arma::repelem(hatInvTauSq.t(), L, 1)+hatInvGammaSq;
		double rateSig = b + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatBeta) % repInvTau));
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
	
	return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
							Rcpp::Named("GS.beta") = gsBeta,
							Rcpp::Named("GS.invTauSq") = gsInvTauSq,
							Rcpp::Named("GS.invGammaSq") = gsInvGammaSq,
							Rcpp::Named("GS.lambda.sq.t") = gsLambdaSqT,
							Rcpp::Named("GS.lambda.sq.g") = gsLambdaSqG,
							Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
							Rcpp::Named("GS.mse") = gsMSE);
}
