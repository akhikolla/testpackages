#include "timeCluster.h"

using namespace Rcpp;
using namespace arma;

RcppExport SEXP calcDistances(SEXP X){

    Rcpp::NumericMatrix XX(X);
    
    int q = XX.ncol();
    int p = XX.nrow();
    mat XXX(XX.begin(), p, q, false);

    mat distances(p,p);
    distances.zeros();
     for(int i=0;i<p;i++){
       for(int j=i;j<p;j++){
      	 for(int k=0;k<q;k++){
	          distances(i,j) = distances(i,j) + (XXX(i,k)-XXX(j,k))*(XXX(i,k)-XXX(j,k)); 
	       }
       }
     }
 
    List res;
    res["distMat"] = sqrt(distances);
    
    return res;
}
