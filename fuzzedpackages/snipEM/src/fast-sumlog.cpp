#include "fast-sumlog.h"
using namespace Rcpp;
using namespace arma;
using namespace std;


/***************************************************/
/*               Function prototypes               */
/***************************************************/
SEXP fast_sumlog(SEXP X, SEXP Lower, SEXP Upper, SEXP N)
{
	try{
		vec x = as<vec>(X);
		double lower = as<double>(Lower);
		double upper = as<double>(Upper);
		int n = as<int>(N);

		uvec index(n); index.fill(1);
		double s = 0; 
		int sind = n; 
		int i = 0; 
		while( i < n-1 && (x(i) < lower || x(i) > upper ) ){
			i++;
		}
		s = x(i); 
		index(i) = 0; 
		sind--; 

		while( sind > 0 ){
			uvec idx = find( index == 1 );
			i = 0;
			while( i < sind - 1 && ( (x(idx(i)) - s) < lower || (x(idx(i)) - s) > upper)){
				i++;
			}
			s += log(1 + exp(x(idx(i)) - s)); 
			index(idx(i)) = 0; 
			sind--; 		
		}
		return wrap(s);	
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	return wrap(NA_REAL);
}
