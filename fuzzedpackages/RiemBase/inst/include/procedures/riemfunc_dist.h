#ifndef RIEMFUNC_DIST_H
#define RIEMFUNC_DIST_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "../include/riemfactory.hpp"

///////////////////////////////////////////////////////////////
inline double riemfunc_dist(arma::mat x, arma::mat y, std::string name){
  if (name=="euclidean"){
    return(euclidean_dist(x,y));
  } else if (name=="sphere"){
    return(sphere_dist(x,y));
  } else if (name=="spd"){
    return(spd_dist(x,y));
  } else if (name=="grassmann"){
    return(grassmann_dist(x,y));
  } else if (name=="stiefel"){
    return(stiefel_dist(x,y));
  } else {
    Rcpp::Rcout << "RiemBase::riemfunc_dist : " << name << " is not yet implemented." << std::endl;
    return(NA_REAL);
  }
}
#endif
