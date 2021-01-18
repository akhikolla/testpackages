#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BRGL (arma::mat xx, arma::vec y, arma::mat W, unsigned int s, unsigned int L, int maxSteps, arma::vec hatAlpha, arma::mat hatBeta, double hatTau, arma::vec hatV, arma::vec hatSg, arma::mat invSigAlpha0, double hatEtaSq, double xi1, double xi2, double r, double a, double b, int progress)
{
	unsigned int n = xx.n_rows, q = W.n_cols;
	arma::mat gsAlpha(maxSteps, q),
			gsBeta(maxSteps, s*L),
			gsV(maxSteps, n),
			gsSg(maxSteps, s);
		
	arma::vec gsEtaSq(maxSteps),
			gsTau(maxSteps),
			gsMSE(maxSteps);
	
	arma::mat varAlpha, tWWoV(q,q), XgXgoV(L,L), temp, varG;
	arma::vec res, RWoV(q), meanAlpha, muV, meanG;
	arma::rowvec RXgoV(L);
	double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), ResSqoV, muS;
	
	std::vector<arma::mat> Xg(s);
	for(unsigned int j=0; j<s; j++){
		Xg[j] = xx.cols((j*L), (j*L+L-1));
	}
	
	for (int k = 0; k < maxSteps; k++) {
		// Rcpp::Rcout << "alpha" << std::endl;
		res = y - xx * arma::vectorise(hatBeta) - xi1*hatV;
		tWWoV = (W.each_col()/hatV).t() * W;
		RWoV = arma::sum(W.each_col()% (res/hatV), 0).t();
		varAlpha = arma::inv(tWWoV*hatTau/xi2Sq + invSigAlpha0);
		meanAlpha = varAlpha * RWoV * hatTau / xi2Sq;
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();
		
		// Rcpp::Rcout << "v" << std::endl;
		res += xi1*hatV;
		lambV = hatTau*xi1Sq/xi2Sq + 2*hatTau;
		muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(res));
		for(unsigned int i = 0; i<n; i++){
			bool flag = true;
			while(flag){
				hatV(i) = 1/rinvGauss(muV(i), lambV);
				if(hatV(i)<=0 || std::isinf(hatV(i)) || std::isnan(hatV(i))){
					if(progress != 0){
						Rcpp::Rcout << "hatV(i) <= 0 or nan or inf" << std::endl;
						Rcpp::checkUserInterrupt();
					}
				}else{
					flag = false;
				}
			}
		}
		res -= xi1*hatV;
		gsV.row(k) = hatV.t();
		
		
		// Rcpp::Rcout << "S" << std::endl;
		for(unsigned int j = 0; j<s; j++){
			muS = std::sqrt(hatEtaSq)/ (arma::norm(hatBeta.col(j)));
			bool flag = true;
			while(flag){
				hatSg(j) = 1/rinvGauss(muS, hatEtaSq);
				if(hatSg(j)<=0 || std::isinf(hatSg(j)) || std::isnan(hatSg(j))){
					if(progress != 0){
						Rcpp::Rcout << "hatSg(j) = " << hatSg(j) << " mu: " << muS << " lamb: " << hatEtaSq << std::endl;
						Rcpp::checkUserInterrupt();
					}
				}else{
					flag = false;
				}
			}
		}
		gsSg.row(k) = hatSg.t();
		
		// res = y - xx * arma::vectorise(hatBeta) - xi1*hatV - W * hatAlpha; //$#@!
		// Rcpp::Rcout << "beta" << std::endl;
		for(unsigned int j=0; j<s; j++){
			res += Xg[j] * hatBeta.col(j);
			XgXgoV = (Xg[j].each_col()/hatV).t() * Xg[j];
			temp = XgXgoV*hatTau/xi2Sq;
			temp.diag() += 1/hatSg(j);
			varG = arma::inv(temp);
			
			RXgoV = arma::sum(Xg[j].each_col()% (res/hatV), 0);
			meanG = varG * RXgoV.t() * hatTau / xi2Sq;
			hatBeta.col(j) = mvrnormCpp(meanG, varG);
			res -= Xg[j] * hatBeta.col(j);
		}
		gsBeta.row(k) = arma::vectorise(hatBeta).t();
		
		// res = y - xx * arma::vectorise(hatBeta) - xi1*hatV - W * hatAlpha; //$#@!
		// Rcpp::Rcout << "tau" << std::endl;
		double shape = a + 3*n/2;
		ResSqoV = arma::accu(arma::square(res)/hatV);
		double rate = b + arma::accu(hatV) + ResSqoV/(2*xi2Sq);
		hatTau = R::rgamma(shape, 1/rate);
		gsTau(k) = hatTau;
		
		
		// Rcpp::Rcout << "eta2Sq" << std::endl;
		double shape2 = (s+s*L)/2+1;
		double rate2 = arma::accu(hatSg)/2 + r;
		hatEtaSq = R::rgamma(shape2, 1/rate2);
		gsEtaSq(k) = hatEtaSq;
		
		
		gsMSE(k) = arma::mean(arma::abs(res));
		if(k % 100 == 0){
			Rcpp::checkUserInterrupt();
		}
		if(progress != 0 && k % progress == 0){
			Rcpp::Rcout << "\nIter." << k << "  MAD: " << gsMSE(k) << std::endl;
		}
	}
	
	return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
							Rcpp::Named("GS.beta") = gsBeta,
							Rcpp::Named("GS.tau") = gsTau,
							Rcpp::Named("GS.v") = gsV,
							Rcpp::Named("GS.s") = gsSg,
							Rcpp::Named("GS.eta2.sq") = gsEtaSq,
							Rcpp::Named("GS.mad") = gsMSE);
}
