
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP file67c49527( SEXP A, SEXP B) ;
}

// definition

SEXP file67c49527( SEXP A, SEXP B ){
BEGIN_RCPP


Rcpp::NumericMatrix FullResp(A);

Rcpp::NumericMatrix Responses(B);



int nex = Responses.nrow();

int nrows = FullResp.nrow();

int item, option, col;





     

         for (int j = 0; j < nrows; j++) {





item = FullResp(j,0);

option = FullResp(j,1);







for(int k = 0; k < nex; k++){





col=k+3;





if(Responses(k,item-1) == option){

FullResp(j,col)=1;

}

else{

FullResp(j,col)=0;

}







         }

     }



     return FullResp;


END_RCPP
}



