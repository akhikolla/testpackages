#include "bbd.h"
using namespace Rcpp;

std::vector<double> prod_mu2_Cpp(const int a, const int A, const int Bp1, const std::vector<double>& mat) {
  
  std::vector<double> res(Bp1*(Bp1+1)/2);
  // res[Trimat(i,j)] = prod(mu2[(i+1):j]), i<j
  // res[Trimax(i,i)] = 1
  
  for (int i = 0; i < (Bp1-1); ++i) {
    res[Trimat(i,i)] = 1;
    res[Trimat(i,i+1)] = mat[(i+1)*(A+1) + a-1];
    for (int j = i+2; j < Bp1; ++j) {
    	res[Trimat(i,j)] = res[Trimat(i,j-1)]*mat[j*(A+1) + a-1];  
      
       // When product of mu2 is infinity, we set it as a large number 
//      if (isinf(res[Trimat(i,j)])) res[Trimat(i,j)] = 1e16;
		}    
  }
  res[Bp1*(Bp1+1)/2 - 1] = 1;
  return(std::move(res));	
}
