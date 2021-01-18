#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"RunLogistic.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;
//using namespace R;

// [[Rcpp::export()]]
arma::mat LogitGrid(arma::mat& x, arma::vec& y, arma::mat& x2, arma::vec& y2, arma::vec lamb1, arma::vec lamb2, arma::vec b, 
                        double r, arma::mat& a, int p, double alpha, char method)
{
	arma::vec btmp(p, fill::zeros);
    arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
	for(unsigned int j=0; j<lamb2.n_elem ; j++){
		for(unsigned int i=0; i<lamb1.n_elem ; i++){
            if(method == 'n'){
                btmp = RunNet(x, y, lamb1(i), lamb2(j), b, r, a, p);
            }else if(method == 'm'){
                btmp = RunMCP(x, y, lamb1(i), b, r, p);
            }else{
                btmp = RunElastic(x, y, lamb1(i), b, alpha, p);
            }
            //btmp = RunNetSurv(xc, xg, y, lamb1(i), lamb2(j), btmp.subvec(0, pc-1), btmp.subvec(pc, pc+p-1), r, a, p, pc, robust);
            CVM(i, j) = validation_logit(x2, y2, btmp);
		}
	}
    return CVM;
}
