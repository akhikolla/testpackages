#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"RunCont.h"
#include"Utilities.h"

using namespace Rcpp;
using namespace arma;
//using namespace R;

// [[Rcpp::export()]]
arma::mat ContGrid(arma::mat& xc, arma::mat& xg, arma::vec& y, arma::mat& x2, arma::vec& y2, arma::vec lamb1, arma::vec lamb2,
					arma::vec bc0, arma::vec bg0, double r, arma::mat& a, int p, int pc, char method)
{
  arma::vec btmp(p, fill::zeros);
  arma::mat CVM(lamb1.n_elem, lamb2.n_elem, fill::zeros);
  
  // std::cout << "xc: "<< xc.n_cols << std::endl;
  // std::cout << "xg: "<< xg.n_cols << std::endl;
  // std::cout << "xc: "<< xc.submat(0, 0, 5, 2) << std::endl;
  // std::cout << "xg: "<< xg.submat(0, 0, 5, 5) << std::endl;
  // std::cout << "y: "<< y.subvec(0, 9) << std::endl;
  // std::cout << "lamb1: "<< lamb1 << std::endl;
  // std::cout << "lamb2: "<< lamb2 << std::endl;
  // std::cout << "p: "<< p << std::endl;
  // std::cout << "pc: "<< pc << std::endl;
  
  for(unsigned int j=0; j<lamb2.n_elem ; j++){
    for(unsigned int i=0; i<lamb1.n_elem ; i++){
	  btmp = RunCont(xc, xg, y, lamb1(i), lamb2(j), bc0, bg0, r, a, p, pc, method);
      CVM(i, j) = validation_LS(x2, y2, btmp);
    }
  }
  return CVM;
}
