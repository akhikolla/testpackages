#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"SurvCD.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;
//using namespace R;
          
// [[Rcpp::export()]]
arma::mat NetGrid(arma::mat& xc, arma::mat& xg, arma::vec& y, arma::mat& x2, arma::vec& y2, arma::vec lamb1, arma::vec lamb2, 
                    arma::vec bc, arma::vec bg, double r, arma::mat& a, int p, int pc, bool robust)
{
	arma::vec bnet((pc+p), fill::zeros);
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	for(unsigned int j=0; j<lamb2.n_elem ; j++){
		for(unsigned int i=0; i<lamb1.n_elem ; i++){
			bnet = RunNetSurv(xc, xg, y, lamb1(i), lamb2(j), bc, bg, r, a, p, pc, robust);
            //bnet = RunNetSurv(xc, xg, y, lamb1(i), lamb2(j), bnet.subvec(0, pc-1), bnet.subvec(pc, pc+p-1), r, a, p, pc, robust);
			if(robust){
				CVM(i, j) = validation_LAD(x2, y2, bnet);
			}else{
				CVM(i, j) = validation_LS(x2, y2, bnet);
			}
		}
        //std::cout << "bc: "<< bnet.subvec(0, pc-1) << std::endl;
	}
	//std::cout << robust << " lamb2: " << lamb2net << " r: " << rnet <<std::endl;
    return CVM;
}


// [[Rcpp::export()]]
arma::vec MCPGrid(arma::mat& xc, arma::mat& xg, arma::vec& y, arma::mat& x2, arma::vec& y2, arma::vec lamb1,
                    arma::vec bc, arma::vec bg, double r, int p, int pc, bool robust)
{
	arma::vec b((pc+p), fill::zeros);
    arma::vec CVM(lamb1.n_elem, fill::zeros);
    for(unsigned int i=0; i<lamb1.n_elem ; i++){
        b =RunMCPSurv(xc, xg, y, lamb1(i), bc, bg, r, p, pc, robust);
        if(robust){
            CVM(i) = validation_LAD(x2, y2, b);
        }else{
            CVM(i) = validation_LS(x2, y2, b);
        }
    }
    return CVM;
}

// [[Rcpp::export()]]
arma::vec LassoGrid(arma::mat& xc, arma::mat& xg, arma::vec& y, arma::mat& x2, arma::vec& y2, arma::vec lamb1,
                    arma::vec bc, arma::vec bg, int p, int pc, bool robust)
{
	arma::vec b((pc+p), fill::zeros);
    arma::vec CVM(lamb1.n_elem, fill::zeros);
    for(unsigned int i=0; i<lamb1.n_elem ; i++){
        b = RunLassoSurv(xc, xg, y, lamb1(i), bc, bg, p, pc, robust);
        if(robust){
            CVM(i) = validation_LAD(x2, y2, b);
        }else{
            CVM(i) = validation_LS(x2, y2, b);
        }
    }
    return CVM;
}
