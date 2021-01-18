#include "local_int_cpp.h"

SEXP local_int_cpp( SEXP posx, SEXP posy, SEXP R1, SEXP R2 ){
    using namespace Rcpp;
    
double R1c = pow(Rcpp::as<double>(R1),2.0);
double R2c = pow(Rcpp::as<double>(R2),2.0);
Rcpp::NumericMatrix xc(posx);
Rcpp::NumericMatrix yc(posy);
int dim = xc.ncol();
int nrow = xc.nrow();
int ncol = yc.nrow();
Rcpp::NumericVector NP(nrow);
std::fill(NP.begin(),NP.end(),0.0);
Rcpp::NumericVector PN(ncol);
std::fill(PN.begin(),PN.end(),0.0);
double dist = 0.0;
    
for (int i = 0; i < nrow; i++) 
{
    for (int j = 0; j < ncol; j++) 
    {
        for(int d = 0; d < dim; d++)
        {
            dist += pow((xc(i,d)-yc(j,d)),2.0);
        }
    if( dist < R1c ) NP[i] += 1;
    if( dist < R2c ) PN[j] += 1;
    dist = 0.0;
    }
}
    
return Rcpp::List::create(Rcpp::Named("x",NP),Rcpp::Named("y",PN));

}
