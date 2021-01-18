//

#include "lmbr.h"



RCPP_MODULE(Clmbr){

	using namespace Rcpp ;


	class_<Clmbr>( "Cpp_Clmbr" )

// expose the constructor
	    .constructor< NumericVector, NumericMatrix, NumericMatrix, int, int, int >()

// expose methods
		.method( "sl3", &Clmbr::sl3R , "SL for th0 by specified method, tolerance" )
		.method( "sl4", &Clmbr::sl4R , "SL for (th0,a0) by specified method, tolerance" )
		.method( "sl5", &Clmbr::sl5R , "SL for th0 by specified method, tolerance and output flag" )
		.method( "sl6", &Clmbr::sl6R , "SL for (th0,a0) by specified method, tolerance and output flag" )
		.method( "ci", &Clmbr::ciR , "confidence interval for theta by specified method" )
		.method( "cr3", &Clmbr::cr3R ,
			"confidence region for (theta,alpha) by specified method and increment" )
		.method( "cr4", &Clmbr::cr4R , 
			"confidence region for (theta,alpha) by specified method and increment, return matrix" )
		.method( "mle", &Clmbr::MLE , "printout maximum likelihood estimates of parameters" )
		.method( "param", &Clmbr::PARAM , "return maximum likelihood estimates of parameters" )
		.method( "sety", &Clmbr::SET_rWy , 
			"reset values for square-root of weights times y-vector" )
	;

}

