#include <Rcpp.h>
#include <functional>
using namespace Rcpp ;

void Rcpp_isomean(
    const NumericMatrix::Row y , 
    double* w, 
    int nn,
    int* k, 
    double* gew, 
    double* ghat
){
  
    int c, j; 
    double neu;

    c = 0;
    k[c] = 0;
    gew[c] = w[0];
    ghat[c] = y[0];
    
    for (j=1; j < nn; j++){  
        c = c+1;
        k[c] = j;
        gew[c] = w[j];
        ghat[c] = y[j];
        
        /* c is at least 1 as nn is > 1 */
        while (ghat[c-1] >= ghat[c]){
            neu = gew[c]+gew[c-1];
            ghat[c-1] = ghat[c-1]+(gew[c]/neu)*(ghat[c]-ghat[c-1]);
            gew[c-1] = neu;
            c = c-1;
            
            if (c==0) break;
        }
    }
    
    while (nn >= 1){
      for (j=k[c]; j < nn; j++){
        ghat[j] = ghat[c];
      }
      nn = k[c];
      c = c-1;
    }
}

// BEGIN_RCPP and END_RCPP bracket the code so that C++ exceptions are passed back to R
extern "C" SEXP misoreg ( SEXP data, SEXP weights ) {
BEGIN_RCPP

    // we don't make modifications to those
    NumericMatrix input(data);
    NumericVector w(weights);
    int n = input.nrow(), m = input.ncol();
    
    NumericMatrix output(n,m);
    
    // since these always have the same size, we just allocate them here once
    IntegerVector k(m) ;
    NumericVector gew(m), ghat(m) ;
    
    for(int i=0; i<n; i++){
        Rcpp_isomean(
            input.row(i), w.begin(), m, 
            k.begin(), gew.begin(), ghat.begin() );
        output.row(i) = ghat ;
    }
    return List::create( 
        _["result"] = output, 
        _["weights"] = w
        ) ;
    
END_RCPP
}

