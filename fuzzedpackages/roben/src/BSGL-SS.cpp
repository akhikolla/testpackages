#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BSGL_SS (arma::mat xx, arma::vec y, arma::mat W, unsigned int s, unsigned int L, int maxSteps, arma::vec hatAlpha, arma::mat hatBg, double hatSigmaSq, arma::mat hatGamma, arma::mat invSigAlpha0, double hatSsq, double hatPi0, double hatPi1, double hatT, double sh0_1, double sh0_0, double sh1_1, double sh1_0, double c, double d, double cutoff, int progress)
{
	unsigned int n = xx.n_rows, q = W.n_cols;
	arma::mat gsAlpha(maxSteps, q),
			gsBeta(maxSteps, s*L),
			gsBg(maxSteps, s*L),
			gsGamma(maxSteps, s*L),
			gsLg(maxSteps, s);
		
	arma::vec gsSsq(maxSteps),
			gsPi0(maxSteps),
			gsPi1(maxSteps),
			gsSigmaSq(maxSteps),
			gsMSE(maxSteps);
	
	arma::mat tWW = W.t()*W, Xg, hatBeta=hatBg%hatGamma, varAlpha, temp, varB;
	arma::vec res, meanAlpha, muV, RXgoV(L), meanB, tBgBg;
	double varG, meanG, lg, t, lgl, tempG;
	
	std::vector<arma::mat> tXgXg(s);
	for(unsigned int j=0; j<s; j++){
		Xg = xx.cols((j*L), (j*L+L-1));
		tXgXg[j] = Xg.t()*Xg;
	}
	
	for (int k = 0; k < maxSteps; k++) {
			
		// Rcpp::Rcout << "alpha" << std::endl;	
		res = y - xx * arma::vectorise(hatBeta);
		varAlpha = arma::inv_sympd(tWW/hatSigmaSq + invSigAlpha0);
		meanAlpha = varAlpha * (W.t() * res/hatSigmaSq);
		hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
		res -= W * hatAlpha;
		gsAlpha.row(k) = hatAlpha.t();
		
		
		// Rcpp::Rcout << "bg" << std::endl;
		double nonzeroBg = 0;
		for(unsigned int j=0; j<s; j++){
			res += xx.cols((j*L), (j*L+L-1)) * hatBeta.col(j);
			
			temp = arma::diagmat(hatGamma.col(j)) * tXgXg[j] * arma::diagmat(hatGamma.col(j)) / hatSigmaSq;
			temp.diag() += 1;
			varB = arma::inv_sympd(temp);
			
			RXgoV = arma::diagmat(hatGamma.col(j)) * xx.cols((j*L), (j*L+L-1)).t() * res / hatSigmaSq;
			meanB = varB * RXgoV;
			
			double lg_temp = arma::as_scalar(exp(-0.5*(RXgoV.t()*varB*RXgoV))*std::sqrt(arma::det(temp)));
			lg = hatPi0/(hatPi0+(1-hatPi0)*lg_temp);
			gsLg(k, j) = lg;
			t = R::runif(0, 1);
			if(t<lg){
				hatBg.col(j) = mvrnormCpp(meanB, varB);
				nonzeroBg++;
			}else{
				hatBg.col(j).zeros();
			}
			hatBeta.col(j) = hatBg.col(j) % hatGamma.col(j);
			res -= xx.cols((j*L), (j*L+L-1)) * hatBeta.col(j);
		}
		gsBg.row(k) = arma::vectorise(hatBg).t();
		
		// Rcpp::Rcout << "Gamma" << std::endl;
		double nonzeroGam = 0;
		for(unsigned int j = 0; j<s; j++){
			for(unsigned int l = 0; l<L; l++){
				res += xx.col(j*L+l) * hatBeta(l,j);
				
				tempG = arma::accu(arma::square(xx.col(j*L+l))) * std::pow(hatBg(l,j),2)/hatSigmaSq + 1/hatSsq;
				varG = 1/tempG;
				double sigG = std::sqrt(varG);
				
				double RXbToV = arma::sum(xx.col(j*L+l) % res) * hatBg(l,j)/hatSigmaSq;
				meanG = varG * RXbToV;
				
				double lgl_temp = 0.5*std::sqrt(hatSsq)*std::exp(-0.5*varG*std::pow(RXbToV,2))/sigG/R::pnorm((meanG/sigG),0,1,1,0);
				// double lgl_temp = 0.5*std::sqrt(hatSsq)*exp(-0.5*varG*pow(RXbToV,2))/sigG/R::pnorm(0, meanG, sigG, 0, 0);
				lgl = hatPi1/(hatPi1+(1-hatPi1)*lgl_temp);
				t = R::runif(0, 1);
				if(t<lgl){
					hatGamma(l,j) = rtnorm0(meanG, sigG);
					nonzeroGam++;
				}else{
					hatGamma(l,j) = 0;
				}
				hatBeta(l,j) = hatBg(l,j) * hatGamma(l,j);
				res -= xx.col(j*L+l) * hatBeta(l,j);
			}
		}
		gsGamma.row(k) = arma::vectorise(hatGamma).t();
		gsBeta.row(k) = arma::vectorise(hatBeta).t();
		
		nonzeroGam = arma::accu(hatGamma > cutoff); //////
		
		// Rcpp::Rcout << "Ssq" << std::endl;
		double shape1 = nonzeroGam/2+1;
		double rate1 = arma::accu(arma::square(hatGamma))/2 + hatT;
		hatSsq = 1/R::rgamma(shape1, 1/rate1);
		gsSsq(k) = hatSsq;
		
		if((k+1) % 10 == 0){
			Rcpp::checkUserInterrupt();
			hatT = 1/arma::mean(1/gsSsq.subvec(k-9,k));
		}
		
		// Rcpp::Rcout << "pi0" << std::endl;
		// tBgBg = arma::sum(arma::square(hatBg), 0).t();
		// double shape0_1 = sh0_1 + arma::accu(tBgBg != 0);
		// double shape0_2 = sh0_0 + arma::accu(tBgBg == 0);
		double shape0_1 = sh0_1 + nonzeroBg;
		double shape0_2 = sh0_0 + (s - nonzeroBg);
		hatPi0 = R::rbeta(shape0_1, shape0_2);
		gsPi0(k) = hatPi0;
		
		// Rcpp::Rcout << "pi1" << std::endl;
		double shape1_1 = sh1_1 + nonzeroGam;
		double shape1_2 = sh1_0 + s*L - nonzeroGam;
		hatPi1 = R::rbeta(shape1_1, shape1_2);
		gsPi1(k) = hatPi1;
		
		// Rcpp::Rcout << "sigma.sq" << std::endl;
		double shapeSig = c + n/2;
		double rateSig = d + 0.5*arma::accu(arma::square(res));
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
							Rcpp::Named("GS.bg") = gsBg,
							Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
							Rcpp::Named("GS.gamma") = gsGamma,
							Rcpp::Named("GS.lg") = gsLg,
							Rcpp::Named("GS.s") = gsSsq,
							Rcpp::Named("GS.pi0") = gsPi0,
							Rcpp::Named("GS.pi1") = gsPi1,
							Rcpp::Named("GS.mse") = gsMSE);
}
