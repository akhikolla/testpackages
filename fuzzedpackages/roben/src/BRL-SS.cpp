#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BRL_SS (arma::mat xx, arma::vec y, arma::mat W, int maxSteps, arma::vec hatAlpha, arma::vec hatBeta, double hatTau, arma::vec hatV, arma::vec hatSg, arma::mat invSigAlpha0, double hatPi, double hatEtaSq, double xi1, double xi2, double r1, double a, double b, double sh1, double sh0, int progress)
{
	unsigned int n = xx.n_rows, s = xx.n_cols, q = W.n_cols;
	arma::mat gsAlpha(maxSteps, q),
			gsBeta(maxSteps, s),
			gsV(maxSteps, n),
			gsSg(maxSteps, s);
		
	arma::vec gsEtaSq(maxSteps),
			gsPi(maxSteps),
			gsTau(maxSteps),
			gsMSE(maxSteps);
	
	arma::mat varAlpha, tWWoV(q,q), temp;
	arma::vec res, RWoV(q), meanAlpha, muV, muS;
	double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), XgXgoV, RXjToV, meanG, varG, lj, t, ResSqoV;
	
	
	for (int k = 0; k < maxSteps; k++) {
		// Rcpp::Rcout << "alpha" << std::endl;
		res = y - xx * hatBeta - xi1*hatV;
		tWWoV = (W.each_col()/hatV).t() * W;
		RWoV = arma::sum(W.each_col()% (res/hatV), 0).t();
		varAlpha = arma::inv_sympd(tWWoV*hatTau/xi2Sq + invSigAlpha0);
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
		muS = std::sqrt(hatEtaSq)/ arma::abs(hatBeta);
		for(unsigned int j = 0; j<s; j++){
			if(hatBeta(j) == 0){
				hatSg(j) = R::rexp(2/hatEtaSq);
			}else{
				bool flag = true;
				while(flag){
					hatSg(j) = 1/rinvGauss(muS(j), hatEtaSq);
					if(hatSg(j)<=0 || std::isinf(hatSg(j)) || std::isnan(hatSg(j))){
						if(progress != 0){
							Rcpp::Rcout << "hatSg(j)： " << hatSg(j) << std::endl; 
							Rcpp::checkUserInterrupt();
						}
					}else{
						flag = false;
					}
				}
			}
			
		}
		gsSg.row(k) = hatSg.t();
		
		// Rcpp::Rcout << "beta" << std::endl;
		for(unsigned int j=0; j<s; j++){
			res += xx.col(j) * hatBeta(j);
			XgXgoV = arma::as_scalar((xx.col(j)/hatV).t() * xx.col(j));
			varG = 1/(XgXgoV*hatTau/xi2Sq + 1/hatSg(j));
			
			RXjToV = arma::sum(xx.col(j) % res / hatV)* hatTau/xi2Sq;
			meanG = varG * RXjToV;
			
			double lj_temp = std::sqrt(hatSg(j))*std::exp(-0.5*varG*pow(RXjToV,2))/std::sqrt(varG);
			lj = hatPi/(hatPi+(1-hatPi)*lj_temp);
			t = R::runif(0, 1);
			if(t<lj){
				hatBeta(j) = R::rnorm(meanG, sqrt(varG));
			}else{
				hatBeta(j) = 0;
			}
			res -= xx.col(j) * hatBeta(j);
		}
		gsBeta.row(k) = hatBeta.t();
		
		
		// Rcpp::Rcout << "etaSq" << std::endl;
		double rate2 = arma::accu(hatSg)/2 + r1;
		hatEtaSq = R::rgamma(s+1, 1/rate2);
		gsEtaSq(k) = hatEtaSq;
		
		// Rcpp::Rcout << "pi" << std::endl;
		double shape1 = sh1 + arma::accu(hatBeta != 0);
		double shape2 = sh0 + arma::accu(hatBeta == 0);
		hatPi = R::rbeta(shape1, shape2);
		gsPi(k) = hatPi;
		
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
							Rcpp::Named("GS.s") = gsSg,
							Rcpp::Named("GS.pi") = gsPi,
							Rcpp::Named("GS.eta2.sq") = gsEtaSq,
							Rcpp::Named("GS.mad") = gsMSE);
}
