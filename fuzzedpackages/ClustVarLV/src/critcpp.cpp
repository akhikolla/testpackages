#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP critcpp(SEXP a, SEXP b){
     NumericMatrix XtX(a);
     NumericVector crit(b);
     int p = XtX.ncol();
          NumericMatrix critstep(p,p);
     NumericMatrix deltamin(p,p);
     List lst(2);
         for (int i = 0; i < (p-1); i++){
            for (int j = i+1; j < p; j++){
               critstep(i,j) = (XtX(i,i)+XtX(j,j)+sqrt((XtX(i,i)+XtX(j,j))*(XtX(i,i)+XtX(j,j))-4*(XtX(i,i)*XtX(j,j)-XtX(i,j)*XtX(j,i))))/(2);
               deltamin(i,j) = crit[i]+crit[j]-critstep(i,j);
             }
          }
          lst[0] = critstep;
          lst[1] = deltamin;
          return lst;
}
