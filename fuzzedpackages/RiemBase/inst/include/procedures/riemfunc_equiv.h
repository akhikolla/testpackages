#ifndef RIEMFUNC_EQUIV_H
#define RIEMFUNC_EQUIV_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "../include/riemfactory.hpp"

///////////////////////////////////////////////////////////////
inline arma::vec riemfunc_equiv(arma::mat x, int m, int n, std::string name){
  if (name=="euclidean"){
    return(euclidean_equiv(x,m,n));
  } else if (name=="sphere"){
    return(sphere_equiv(x,m,n));
  } else if (name=="spd"){
    return(spd_equiv(x,m,n));
  } else if (name=="grassmann"){
    return(grassmann_equiv(x,m,n));
  } else if (name=="stiefel"){
    return(stiefel_equiv(x,m,n));
  } else {
    Rcpp::Rcout << "RiemBase::riemfunc_equiv : " << name << " is not yet implemented." << std::endl;
    Rcpp::stop(""); 
  }
}
#endif
