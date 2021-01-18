
#include <Rcpp.h>
#include "jmc.h"
#include "jmo.h"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List  jmc_main(SEXP k, SEXP n1,SEXP p1,SEXP p2, SEXP maxl, SEXP p1a, SEXP maxiterations, SEXP point,SEXP xs,SEXP ws,SEXP yfile, SEXP cfile, SEXP mfile,SEXP trace)
  {
	Rcpp::List result;  
  try {
     
	 result=jmcspace::jmc_cmain( as<int> (k), as<int> (n1), as<int> (p1), as<int> (p2), 
           as<int> (maxl), as<int> (p1a), as<int>(maxiterations),as<int>(point), as<std::vector<double> >(xs),as<std::vector<double> >(ws), as<std::string> (yfile), 
           as<std::string> (cfile),as<std::string>(mfile),as<int> (trace));
   if(Rf_isNull(result)){
	    throw std::range_error("Possible files reading or format errors");
        }
     return result;			   
	} catch(std::exception &ex) {	
	forward_exception_to_r(ex);
    } catch(...) { 
	::Rf_error("c++ exception (unknown reason)"); 
    }
    return R_NilValue;             // not reached
  
}

// [[Rcpp::export]]
Rcpp::List  jmo_main(SEXP k, SEXP n1,SEXP p1,SEXP p2, SEXP p1a, SEXP bq,SEXP K_num, SEXP j_max,  SEXP point,SEXP xs,SEXP ws,SEXP betas,SEXP thetas,SEXP maxiterations, SEXP yfile, SEXP cfile, SEXP mfile,SEXP trace)
  {
	Rcpp::List result;  
  try {
     
	 result=jmospace::jmo_cmain( as<int> (k), as<int> (n1), as<int> (p1), as<int> (p2), 
           as<int> (p1a),  as<int> (bq), as<int> (K_num), as<int> (j_max),as<int>(point),as<std::vector<double> >(xs),as<std::vector<double> >(ws), as<std::vector<double> >(betas),as<std::vector<double> >(thetas),
		   as<int>(maxiterations), as<std::string> (yfile), 
           as<std::string> (cfile),as<std::string>(mfile),as<int> (trace));
   if(Rf_isNull(result)){
	    throw std::range_error("Possible files reading or format errors");
        }
     return result;			   
	} catch(std::exception &ex) {	
	forward_exception_to_r(ex);
    } catch(...) { 
	::Rf_error("c++ exception (unknown reason)"); 
    }
    return R_NilValue;             // not reached
  
}






