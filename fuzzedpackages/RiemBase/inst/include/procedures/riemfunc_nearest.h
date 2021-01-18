#ifndef RIEMFUNC_NEAREST_H
#define RIEMFUNC_NEAREST_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "../include/riemfactory.hpp"

///////////////////////////////////////////////////////////////
inline arma::mat riemfunc_nearest(arma::mat x, std::string name){
  if (name=="euclidean"){
    return(euclidean_nearest(x));
  } else if (name=="sphere"){
    return(sphere_nearest(x));
  } else if (name=="spd"){
    return(spd_nearest(x));
  } else if (name=="stiefel"){
    return(stiefel_nearest(x));
  } else if (name=="grassmann"){
    return(grassmann_nearest(x));
  } else {
    Rcpp::Rcout << "RiemBase::riemfunc_nearest : " << name << " is not yet implemented." << std::endl;
    Rcpp::stop("");
  }
}

#endif
