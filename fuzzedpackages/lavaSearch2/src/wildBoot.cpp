// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat wildBoot_cpp (const arma::mat& iid, std::vector< arma::uvec > lsIndexModel,
						int nSample, int nObs, int nModel, int p) {
  // for modelsearch2 (sampleMaxDistChi2)
  arma::colvec weight(nObs);
  arma::mat wiid(nSample,p);
  arma::rowvec iScore(p);
  arma::mat iStat(nSample,nModel);
  
  for(int iSample=0; iSample<nSample; iSample++){
	wiid = iid;
	weight.randn(); // generate weights under the null
	wiid.each_col() %= weight;
	iScore = pow(sum(wiid,0),2);
	
	for(int iModel=0; iModel<nModel; iModel++){
	  arma::rowvec tempo = iScore.cols(lsIndexModel[iModel]);
	  iStat(iSample,iModel) = sum(tempo);
	}
		
  }

  return(iStat);
}                  

