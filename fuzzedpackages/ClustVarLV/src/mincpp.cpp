#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP mincpp(SEXP a){
          NumericMatrix X(a);
          int p = X.ncol();
          List lst(3);
          double minv = 1000000000;
          int mini = 0;
          int minj = 1;
     for (int i = 0; i < (p-1); i++){
          for (int j = i+1; j < p; j++){
              if(X(i,j) < minv){
                      minv=X(i,j);
                      mini = i+1;
                      minj = j+1;
               }
         }
     }
     lst[0] = mini;
     lst[1] = minj;
     lst[2] = minv;
     return lst;
   }
