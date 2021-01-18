#ifndef RIEMFUNC_EXTDIST_H
#define RIEMFUNC_EXTDIST_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "../include/riemfactory.hpp"

///////////////////////////////////////////////////////////////
inline double riemfunc_extdist(arma::mat x, arma::mat y, std::string name){
  if (name=="euclidean"){
    return(euclidean_extdist(x,y));
  } else if (name=="sphere"){
    return(sphere_extdist(x,y));
  } else if (name=="spd"){
    return(spd_extdist(x,y));
  } else if (name=="grassmann"){
    return(grassmann_extdist(x,y));
  } else if (name=="stiefel"){
    return(stiefel_extdist(x,y));
  } else {
    Rcpp::Rcout << "RiemBase::riemfunc_extdist : " << name << " is not yet implemented." << std::endl;
    return(NA_REAL);
  }
}
#endif
