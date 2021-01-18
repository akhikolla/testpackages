#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BRSGL_SS (arma::mat xx, arma::vec y, arma::mat W, unsigned int s, unsigned int L, int maxSteps, arma::vec hatAlpha, arma::mat hatBg, double hatTau, arma::vec hatV, arma::mat hatGamma, arma::mat invSigAlpha0, double hatSsq, double hatPi0, double hatPi1, double xi1, double xi2, double hatT, double a, double b, double sh0_1, double sh0_0, double sh1_1, double sh1_0, double cutoff, int progress)
{
	unsigned int n = xx.n_rows, q = W.n_cols;
	arma::mat gsAlpha(maxSteps, q),
			gsBeta(maxSteps, s*L),
			gsBg(maxSteps, s*L),
			gsV(maxSteps, n),
			gsGamma(maxSteps, s*L),
			gsLg(maxSteps, s);
		
	arma::vec gsSsq(maxSteps),
			gsPi0(maxSteps),
			gsPi1(maxSteps),
			gsTau(maxSteps),
			gsMSE(maxSteps);
	
	arma::mat hatBeta=hatBg%hatGamma, varAlpha, tWWoV(q,q), XgXgoV(L,L), temp, varB;
	arma::vec res, RWoV(q), meanAlpha, muV, RXgoV(L), meanB, tBgBg;
	double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), varG, meanG, lg, t, lgl, ResSqoV, tempG;
	
	
	for (int k = 0; k < maxSteps; k++) {
		// Rcpp::Rcout << k << std::endl;
		// Rcpp::Rcout << "alpha" << std::endl;
		res = y - xx * arma::vectorise(hatBeta) - xi1*hatV;		
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
		muV = std::sqrt(xi1Sq+2*xi2Sq) / arma::abs(res);
		// if(muV.has_nan()){
			// std::string error = std::string("v: NaN in nu of rinvGauss");
			// throw std::runtime_error(error);
		// }
		for(unsigned int i = 0; i<n; i++){
			bool flag = true;
			while(flag){
				hatV(i) = 1/rinvGauss(muV(i), lambV);
				if(hatV(i)<=0 || std::isinf(hatV(i)) || std::isnan(hatV(i))){
					if(progress != 0){
						Rcpp::Rcout << "hatV(i) = " << hatV(i) << " muV(i): " << muV(i) << " (res)(i): " << res(i) << " lambV: " << lambV << std::endl;
						Rcpp::checkUserInterrupt();
					}
				}else{
					flag = false;
				}
			}
		}
		res -= xi1*hatV;
		gsV.row(k) = hatV.t();
		
		// Rcpp::Rcout << "bg" << std::endl; // group level
		double nonzeroBg = 0;
		for(unsigned int j=0; j<s; j++){
			res += xx.cols((j*L), (j*L+L-1)) * hatBeta.col(j);
			
			XgXgoV = (xx.cols((j*L), (j*L+L-1)).each_col()/hatV).t() * xx.cols((j*L), (j*L+L-1));
			temp = arma::diagmat(hatGamma.col(j)) * XgXgoV*hatTau/xi2Sq * arma::diagmat(hatGamma.col(j));
			temp.diag() += 1;
			// varB = arma::inv_sympd(temp);
			// if(!temp.is_symmetric()){
				// Rcpp::Rcout << "temp is not symmetric:\n" << std::endl;
			// }
			varB = arma::inv_sympd(temp);
			RXgoV = arma::sum(xx.cols((j*L), (j*L+L-1)).each_col()% (res/hatV), 0).t() % hatGamma.col(j);
			RXgoV *=  hatTau/xi2Sq;
			meanB = varB * RXgoV;
			
			double lg_temp = arma::as_scalar(exp(-0.5*(RXgoV.t()*varB*RXgoV)))*std::sqrt(arma::det(temp));
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
		
		// Rcpp::Rcout << "Gamma" << std::endl; // individual level
		double nonzeroGam = 0;
		for(unsigned int j = 0; j<s; j++){
			for(unsigned int l = 0; l<L; l++){
				res += xx.col(j*L+l) * hatBeta(l,j);
				
				tempG = arma::accu(arma::square(xx.col(j*L+l))/hatV) * std::pow(hatBg(l,j),2) * hatTau/xi2Sq + 1/hatSsq;
				varG = 1/tempG;
				double sigG = std::sqrt(varG);
				
				double RXbToV = arma::sum(xx.col(j*L+l) % (res/hatV)) * hatBg(l,j) * hatTau/xi2Sq;
				meanG = varG * RXbToV;
				
				double lgl_temp = 0.5*std::sqrt(hatSsq)*std::exp(-0.5*varG*pow(RXbToV,2))/sigG/R::pnorm((meanG/sigG),0,1,1,0);
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
		// double shape1 = arma::accu(hatGamma != 0)/2+1;
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
		
		// Rcpp::Rcout << "tau" << std::endl;
		double shape = a + 3*n/2;
		ResSqoV = arma::accu(arma::square(res)/hatV);
		double rate = b + arma::accu(hatV) + ResSqoV/(2*xi2Sq);
		hatTau = R::rgamma(shape, 1/rate);
		gsTau(k) = hatTau;
		
		gsMSE(k) = arma::mean(arma::abs(res));
		if(progress != 0 && k % progress == 0){
			Rcpp::Rcout << "\nIteration: " << k << std::endl;
			Rcpp::Rcout << "  MAD      : " << gsMSE(k) << std::endl;
			Rcpp::Rcout << "  shape0_1 : " << shape0_1 << "  shape0_0 : " << shape0_2 << std::endl;
			Rcpp::Rcout << "  shape1_1 : " << shape1_1 << "  shape1_0 : " << shape1_2 << std::endl;
		}
	}
	
	return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
							Rcpp::Named("GS.beta") = gsBeta,
							Rcpp::Named("GS.bg") = gsBg,
							Rcpp::Named("GS.tau") = gsTau,
							Rcpp::Named("GS.v") = gsV,
							Rcpp::Named("GS.gamma") = gsGamma,
							Rcpp::Named("GS.lg") = gsLg,
							Rcpp::Named("GS.s") = gsSsq,
							Rcpp::Named("GS.pi0") = gsPi0,
							Rcpp::Named("GS.pi1") = gsPi1,
							Rcpp::Named("GS.mad") = gsMSE);
}
