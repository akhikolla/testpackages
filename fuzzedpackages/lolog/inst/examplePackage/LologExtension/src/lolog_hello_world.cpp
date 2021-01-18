// [[Rcpp::depends(lolog)]]
#include "lolog_hello_world.h"
#include <lolog.h>

//' An basic example of a function in C++ using lolog
//' @return a list of a character vector, a numeric vector, and an lolog DirectedNet
//' @examples
//' lolog_hello_world()
//'
//' #MinDegree is an new statistic defined by this package
//' if(require(network)){
//' 	data(flo)
//' 	nflo<-as.BinaryNet(network(flo,directed=FALSE) )
//' 	fit <- lolog(nflo ~ edges() + minDegree(3),verbose=0)
//' 	summary(fit)
//' }
// [[Rcpp::export]]
Rcpp::List lolog_hello_world(){
    using namespace Rcpp;
    
    IntegerMatrix tmp(0,2);
    lolog::DirectedNet net(tmp,20); 
    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y, wrap(net) ) ;
    
    return z ;
}
