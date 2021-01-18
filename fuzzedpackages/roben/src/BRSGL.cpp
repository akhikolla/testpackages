#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BRSGL (arma::mat xx, arma::vec y, arma::mat W, unsigned int s, unsigned int L, int maxSteps, arma::vec hatAlpha, arma::mat hatBeta, double hatTau, arma::vec hatV, arma::vec hatSg, arma::mat hatGamma, arma::mat invSigAlpha0, double hatEta1Sq, double hatEta2Sq, double xi1, double xi2, double s1, double s2, double r1, double r2, double a, double b, int progress)
{
	unsigned int n = xx.n_rows, q = W.n_cols;
	arma::mat gsAlpha(maxSteps, q),
			gsBeta(maxSteps, s*L),
			gsV(maxSteps, n),
			gsGamma(maxSteps, s*L),
			gsSg(maxSteps, s);
		
	arma::vec gsEta1Sq(maxSteps),
			gsEta2Sq(maxSteps),
			gsTau(maxSteps),
			gsMSE(maxSteps);
	
	arma::mat varAlpha, tWWoV(q,q), XgXgoV(L,L), temp, varG;
	arma::vec res, RWoV(q), meanAlpha, muV, meanG, muG;
	arma::rowvec RXgoV(L);
	double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), ResSqoV, muS;
	
	std::vector<arma::mat> Xg(s);
	for(unsigned int j=0; j<s; j++){
		Xg[j] = xx.cols((j*L), (j*L+L-1));
	}
	
	std::vector<arma::mat> tWW(n);
	for(unsigned int i=0; i<n; i++){
		tWW[i] = W.row(i).t()*W.row(i);
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
					if(progress != 0) Rcpp::Rcout << "hatV(i) <= 0 or nan or inf" << std::endl; 
					Rcpp::checkUserInterrupt();
				}else{
					flag = false;
				}
			}
		}
		res -= xi1*hatV;
		gsV.row(k) = hatV.t();
		
		// Rcpp::Rcout << "beta" << std::endl;
		for(unsigned int j=0; j<s; j++){
			res += Xg[j] * hatBeta.col(j);
			XgXgoV = (Xg[j].each_col()/hatV).t() * Xg[j];
			RXgoV = arma::sum(Xg[j].each_col()% (res/hatV), 0);
			
			temp = XgXgoV*hatTau/xi2Sq;
			temp.diag() += 1/hatGamma.col(j) + 1/hatSg(j);
			varG = arma::inv(temp);
			meanG = varG * hatTau * RXgoV.t() / xi2Sq;
			
			hatBeta.col(j) = mvrnormCpp(meanG, varG);
			res -= Xg[j] * hatBeta.col(j);
		}
		gsBeta.row(k) = arma::vectorise(hatBeta).t();
		
		// Rcpp::Rcout << "Gamma" << std::endl;
		for(unsigned int j = 0; j<s; j++){
			muG = std::sqrt(hatEta1Sq)/ arma::abs(hatBeta.col(j));
			for(unsigned int l = 0; l<L; l++){
				bool flag = true;
				while(flag){
					hatGamma(l,j) = 1/rinvgaussian(muG(l), hatEta1Sq);
					if(hatGamma(l,j)<=0 || std::isinf(hatGamma(l,j)) || std::isnan(hatGamma(l,j))){
						if(progress != 0){
							Rcpp::Rcout << "hatGamma(l,j)： " << hatGamma(l,j) << std::endl; 
							Rcpp::checkUserInterrupt();
						}
					}else{
						flag = false;
					}
				}
			}
		}
		gsGamma.row(k) = arma::vectorise(hatGamma).t();
		
		// Rcpp::Rcout << "S" << std::endl;
		for(unsigned int j = 0; j<s; j++){
			muS = std::sqrt(hatEta2Sq)/ (arma::norm(hatBeta.col(j)));
			hatSg(j) = 1/rinvgaussian(muS, hatEta2Sq);
		}
		gsSg.row(k) = hatSg.t();
		
		// Rcpp::Rcout << "eta1Sq" << std::endl;
		double shape1 = s*L+s1;
		double rate1 = arma::accu(hatGamma)/2 + r1;
		hatEta1Sq = R::rgamma(shape1, 1/rate1);
		gsEta1Sq(k) = hatEta1Sq;
		
		// Rcpp::Rcout << "eta2Sq" << std::endl;
		double shape2 = s/2+s2;
		double rate2 = arma::accu(hatSg)/2 + r2;
		hatEta2Sq = R::rgamma(shape2, 1/rate2);
		gsEta2Sq(k) = hatEta2Sq;
		
		// Rcpp::Rcout << "tau" << std::endl;
		double shape = a + 3*n/2;
		ResSqoV = arma::accu(arma::square(res)/hatV);
		double rate = b + arma::accu(hatV) + ResSqoV/(2*xi2Sq);
		hatTau = R::rgamma(shape, 1/rate);
		gsTau(k) = hatTau;

		
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
							Rcpp::Named("GS.gamma") = gsGamma,
							Rcpp::Named("GS.s") = gsSg,
							Rcpp::Named("GS.eta1.sq") = gsEta1Sq,
							Rcpp::Named("GS.eta2.sq") = gsEta2Sq,
							Rcpp::Named("GS.mad") = gsMSE);
}
