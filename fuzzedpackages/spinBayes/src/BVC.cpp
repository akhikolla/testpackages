#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BVCStr (arma::mat xx, arma::vec y, arma::mat W, arma::mat Ex, unsigned int s, unsigned int q, int maxSteps, arma::vec hatM, arma::vec hatR0, arma::vec hatRStar, arma::vec hatAlpha, arma::vec hatZeta, arma::vec hatInvSigM0, arma::vec hatInvTauSq0, arma::vec hatInvTauSqStar, arma::mat invSigAlpha0, arma::vec hatInvTauSqZeta, double hatLambdaSq0, double hatLambdaSqStar, double hatLambdaSqZeta, double hatSigmaSq, double a0, double b0, double aStar, double bStar, double alpha, double gamma, int progress)
{
	unsigned int L = q-1, n = xx.n_rows, clc = W.n_cols;
	arma::mat gsM(maxSteps, q),
			gsAlpha(maxSteps, clc),
			gsR0(maxSteps, s),
			gsZeta(maxSteps, s),
			gsRStar(maxSteps, s*(q-1)),
			gsInvTauSq0(maxSteps, s),
			gsInvTauSqZeta(maxSteps, s),
			gsInvTauSqStar(maxSteps, s);
		
	arma::vec gsLambda0(maxSteps),
			gsLambdaZeta(maxSteps),
			gsLambdaStar(maxSteps),
			gsSigmaSq(maxSteps);

	arma::mat Bm = xx.cols(0, q-1), B0 = xx.cols(q, q+s-1), Br = xx.cols(q+s, xx.n_cols-1);
	arma::mat tBmBm = Bm.t()*Bm, tB0B0 = B0.t()*B0, tWW = W.t()*W, tExEx = Ex.t()*Ex;
	arma::vec tB0B0Diag = tB0B0.diag();
	arma::vec tExExDiag = tExEx.diag();
	arma::mat invSigM0 = arma::diagmat(hatInvSigM0);
	
	arma::mat Xr, varM, varAlpha, varRs, tempS, matRStar;
	arma::vec res, BrjtRes, meanM,  meanAlpha, meanRs, tRsRs, repInvTau, muInvTauSq0, muInvTauSqStar, muInvTauSqZeta; // mu_m, mu_alpha, 
	double temp0, meanR0, varR0, tempZ, meanZ, varZ, B0jtRes, lInvTauSq0, ExjtRes, lInvTauSqStar, lInvTauSqZeta;
	
	std::vector<arma::mat> tBrBr(s);
	for(unsigned int j=0; j<s; j++){
		Xr = Br.cols((j*L), (j*L+L-1));
		tBrBr[j] = Xr.t()*Xr;
	}
	
	for (int k = 0; k < maxSteps; k++) {		
		// m|y, r0, r.star
		varM = arma::inv(tBmBm/hatSigmaSq + invSigM0);
		res = y - (W * hatAlpha + Ex * hatZeta + B0 * hatR0 + Br * hatRStar);
		meanM = varM * (Bm.t() * res/hatSigmaSq);
		hatM = mvrnormCpp(meanM, varM);
		res -= Bm * hatM;
		gsM.row(k) = hatM.t();
		
		// alpha|y, m, r0, r.star
		varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
		res += W * hatAlpha;
		meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();
		
		for(unsigned int j=0; j<s; j++){
			temp0 = 1/(tB0B0Diag(j) + hatInvTauSq0(j));
			varR0 = hatSigmaSq * temp0;
			res += B0.col(j) * hatR0(j);
			B0jtRes = arma::as_scalar(B0.col(j).t() * res);
			meanR0 = temp0 * B0jtRes;
			hatR0(j) = R::rnorm(meanR0, sqrt(varR0));
			res -= B0.col(j) * hatR0(j);
			
			tempZ = 1/(tExExDiag(j) + hatInvTauSqZeta(j));
			varZ = hatSigmaSq * tempZ;
			res += Ex.col(j) * hatZeta(j);
			ExjtRes = arma::as_scalar(Ex.col(j).t() * res);
			meanZ = tempZ * ExjtRes;
			hatZeta(j) = R::rnorm(meanZ, sqrt(varZ));
			res -= Ex.col(j) * hatZeta(j);
			
			tempS = tBrBr[j];
			tempS.diag() += hatInvTauSqStar(j);
			tempS = arma::inv(tempS);
			varRs = hatSigmaSq * tempS;
			res += Br.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
			BrjtRes = Br.cols((j*L), (j*L+L-1)).t() * res;
			meanRs = tempS * BrjtRes;
			hatRStar.subvec((j*L), (j*L+L-1)) = mvrnormCpp(meanRs, varRs);
			res -= Br.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
		}
		gsR0.row(k) = hatR0.t();
		gsZeta.row(k) = hatZeta.t();
		gsRStar.row(k) = hatRStar.t();
		
		// sigma.sq|
		double shapeSig = alpha + (n+2*s+s*L)/2;
		repInvTau = arma::vectorise(arma::repelem(hatInvTauSqStar.t(), L, 1), 0);
		double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatR0) % hatInvTauSq0) + 
									arma::accu(square(hatZeta) % hatInvTauSqZeta) + 
									arma::accu(square(hatRStar) % repInvTau));
		hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
		gsSigmaSq(k) = hatSigmaSq;
		
		// invTAUsq.0|lambda, r0
		lInvTauSq0 = hatLambdaSq0;
		muInvTauSq0 = sqrt(hatLambdaSq0 * hatSigmaSq / square(hatR0));
		for(unsigned int j = 0; j < s; j++){
			hatInvTauSq0(j) = rinvgaussian(muInvTauSq0(j), lInvTauSq0);
		}
		gsInvTauSq0.row(k) = hatInvTauSq0.t();
		
		// invTAUsq.e|
		lInvTauSqZeta = hatLambdaSqZeta;
		muInvTauSqZeta = sqrt(hatLambdaSqZeta * hatSigmaSq / square(hatZeta));
		for(unsigned int j = 0; j < s; j++){
			hatInvTauSqZeta(j) = rinvgaussian(muInvTauSqZeta(j), lInvTauSqZeta);
		}
		gsInvTauSqZeta.row(k) = hatInvTauSqZeta.t();
		
		// invTAUsq.star|lambda.star, r.star
		lInvTauSqStar = L * hatLambdaSqStar;
		matRStar = arma::reshape(hatRStar, L, s);
		tRsRs = sum(square(matRStar), 0).t();
		muInvTauSqStar = sqrt(L * hatLambdaSqStar * hatSigmaSq / tRsRs);
		for(unsigned int j = 0; j<s; j++){
			hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
		}
		gsInvTauSqStar.row(k) = hatInvTauSqStar.t();
		
		// lambda0|invTAUsq.0
		double shape = a0 + s;
		double rate = b0 + arma::accu(1/hatInvTauSq0)/2;
		hatLambdaSq0 = R::rgamma(shape, 1/rate);
		gsLambda0(k) = hatLambdaSq0;
		
		// lambda.zeta|
		double shapeZ = a0 + s;
		double rateZ = b0 + arma::accu(1/hatInvTauSqZeta)/2;
		hatLambdaSqZeta = R::rgamma(shapeZ, 1/rateZ);
		gsLambdaZeta(k) = hatLambdaSqZeta;
		
		// lambda.star|invTAUsq.star
		double shapeS = aStar + s*(L+1)/2;
		double rateS = bStar + L*arma::accu(1/hatInvTauSqStar)/2;
		hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
		gsLambdaStar(k) = hatLambdaSqStar;
		
		if(progress != 0 && k % progress == 0){
			Rcpp::Rcout << "Iteration: " << k << std::endl;
			Rcpp::Rcout << "  mse    : " << arma::accu(arma::square(res))/n << std::endl;
			Rcpp::Rcout << "  sigmaSq: " << hatSigmaSq << std::endl;
		}
	}
	
	return Rcpp::List::create(Rcpp::Named("GS.m") = gsM,
							Rcpp::Named("GS.alpha") = gsAlpha,
							Rcpp::Named("GS.r0") = gsR0,
							Rcpp::Named("GS.zeta") = gsZeta,
							Rcpp::Named("GS.rs") = gsRStar,
							Rcpp::Named("GS.invTAUsq.0") = gsInvTauSq0,
							Rcpp::Named("GS.invTAUsq.zeta") = gsInvTauSqZeta,
							Rcpp::Named("GS.invTAUsq.star") = gsInvTauSqStar,
							Rcpp::Named("GS.lambda.sq.0") = gsLambda0,
							Rcpp::Named("GS.lambda.sq.zeta") = gsLambdaZeta,
							Rcpp::Named("GS.lambda.sq.star") = gsLambdaStar,
							Rcpp::Named("GS.sigma.sq") = gsSigmaSq);
}



// [[Rcpp::export()]]
Rcpp::List BVCStr_NoE (arma::mat xx, arma::vec y, arma::mat W, bool CLIN, unsigned int s, unsigned int q, int maxSteps, arma::vec hatM, arma::vec hatR0, arma::vec hatRStar, arma::vec hatAlpha, arma::vec hatInvSigM0, arma::vec hatInvTauSq0, arma::vec hatInvTauSqStar, arma::mat invSigAlpha0, double hatLambdaSq0, double hatLambdaSqStar, double hatSigmaSq, double a0, double b0, double aStar, double bStar, double alpha, double gamma, int progress)
{
	unsigned int L = q-1, n = xx.n_rows, clc = W.n_cols;
	arma::mat gsM(maxSteps, q),
			gsAlpha(maxSteps, clc),
			gsR0(maxSteps, s),
			gsRStar(maxSteps, s*(q-1)),
			gsInvTauSq0(maxSteps, s),
			gsInvTauSqStar(maxSteps, s);
		
	arma::vec gsLambda0(maxSteps),
			gsLambdaStar(maxSteps),
			gsSigmaSq(maxSteps);

	arma::mat Bm = xx.cols(0, q-1), B0 = xx.cols(q, q+s-1), Br = xx.cols(q+s, xx.n_cols-1);
	arma::mat tBmBm = Bm.t()*Bm, tB0B0 = B0.t()*B0, tWW = W.t()*W; // tExEx = Ex.t()*Ex;
	arma::vec tB0B0Diag = tB0B0.diag();
	arma::mat invSigM0 = arma::diagmat(hatInvSigM0);
	
	arma::mat Xr, varM, varAlpha, varRs, tempS, matRStar;
	arma::vec res, BrjtRes, meanM,  meanAlpha, meanRs, tRsRs, repInvTau, muInvTauSq0, muInvTauSqStar; // mu_m, mu_alpha, 
	double temp0, meanR0, varR0, B0jtRes, lInvTauSq0, lInvTauSqStar;
	
	std::vector<arma::mat> tBrBr(s);
	for(unsigned int j=0; j<s; j++){
		Xr = Br.cols((j*L), (j*L+L-1));
		tBrBr[j] = Xr.t()*Xr;
	}
	
	for (int k = 0; k < maxSteps; k++) {		
		// m|y, r0, r.star
		varM = arma::inv(tBmBm/hatSigmaSq + invSigM0);
		res = y - (W * hatAlpha + B0 * hatR0 + Br * hatRStar);
		meanM = varM * (Bm.t() * res/hatSigmaSq);
		hatM = mvrnormCpp(meanM, varM);
		res -= Bm * hatM;
		gsM.row(k) = hatM.t();
		
		// alpha|y, m, r0, r.star
		if(CLIN){
			varAlpha = arma::inv(tWW/hatSigmaSq + invSigAlpha0);
			res += W * hatAlpha;
			meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
			hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
			res -= W * hatAlpha;
			gsAlpha.row(k) = hatAlpha.t();
		}
		
		for(unsigned int j=0; j<s; j++){
			temp0 = 1/(tB0B0Diag(j) + hatInvTauSq0(j));
			varR0 = hatSigmaSq * temp0;
			res += B0.col(j) * hatR0(j);
			B0jtRes = arma::as_scalar(B0.col(j).t() * res);
			meanR0 = temp0 * B0jtRes;
			hatR0(j) = R::rnorm(meanR0, sqrt(varR0));
			res -= B0.col(j) * hatR0(j);
			
			tempS = tBrBr[j];
			tempS.diag() += hatInvTauSqStar(j);
			tempS = arma::inv(tempS);
			varRs = hatSigmaSq * tempS;
			res += Br.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
			BrjtRes = Br.cols((j*L), (j*L+L-1)).t() * res;
			meanRs = tempS * BrjtRes;
			hatRStar.subvec((j*L), (j*L+L-1)) = mvrnormCpp(meanRs, varRs);
			res -= Br.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
		}
		gsR0.row(k) = hatR0.t();
		gsRStar.row(k) = hatRStar.t();
		
		// sigma.sq|
		double shapeSig = alpha + (n+s+s*L)/2;
		repInvTau = arma::vectorise(arma::repelem(hatInvTauSqStar.t(), L, 1), 0);
		double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatR0) % hatInvTauSq0) + 
									arma::accu(square(hatRStar) % repInvTau));
		hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
		gsSigmaSq(k) = hatSigmaSq;
		
		// invTAUsq.0|lambda, r0
		lInvTauSq0 = hatLambdaSq0;
		muInvTauSq0 = sqrt(hatLambdaSq0 * hatSigmaSq / square(hatR0));
		for(unsigned int j = 0; j < s; j++){
			hatInvTauSq0(j) = rinvgaussian(muInvTauSq0(j), lInvTauSq0);
		}
		gsInvTauSq0.row(k) = hatInvTauSq0.t();
		
		// invTAUsq.star|lambda.star, r.star
		lInvTauSqStar = L * hatLambdaSqStar;
		matRStar = arma::reshape(hatRStar, L, s);
		tRsRs = sum(square(matRStar), 0).t();
		muInvTauSqStar = sqrt(L * hatLambdaSqStar * hatSigmaSq / tRsRs);
		for(unsigned int j = 0; j<s; j++){
			hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
		}
		gsInvTauSqStar.row(k) = hatInvTauSqStar.t();
		
		// lambda0|invTAUsq.0
		double shape = a0 + s;
		double rate = b0 + arma::accu(1/hatInvTauSq0)/2;
		hatLambdaSq0 = R::rgamma(shape, 1/rate);
		gsLambda0(k) = hatLambdaSq0;
		
		// lambda.star|invTAUsq.star
		double shapeS = aStar + s*(L+1)/2;
		double rateS = bStar + L*arma::accu(1/hatInvTauSqStar)/2;
		hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
		gsLambdaStar(k) = hatLambdaSqStar;
	}
	
	return Rcpp::List::create(Rcpp::Named("GS.m") = gsM,
							Rcpp::Named("GS.alpha") = gsAlpha,
							Rcpp::Named("GS.r0") = gsR0,
							Rcpp::Named("GS.rs") = gsRStar,
							Rcpp::Named("GS.invTAUsq.0") = gsInvTauSq0,
							Rcpp::Named("GS.invTAUsq.star") = gsInvTauSqStar,
							Rcpp::Named("GS.lambda.sq.0") = gsLambda0,
							Rcpp::Named("GS.lambda.sq.star") = gsLambdaStar,
							Rcpp::Named("GS.sigma.sq") = gsSigmaSq);
}
