#ifndef RIEMFUNC_LOG_H
#define RIEMFUNC_LOG_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "../include/riemfactory.hpp"

///////////////////////////////////////////////////////////////
inline arma::mat riemfunc_log(arma::mat x, arma::mat y, std::string name){
  if (name=="euclidean"){
    return(euclidean_log(x,y));
  } else if (name=="sphere"){
    return(sphere_log(x,y));
  } else if (name=="spd"){
    return(spd_log(x,y));
  } else if (name=="grassmann"){
    return(grassmann_log(x,y));
  } else if (name=="stiefel"){
    return(stiefel_log(x,y));
  }else {
    Rcpp::Rcout << "RiemBase::riemfunc_log : " << name <<" is not yet implemented." << std::endl;
    Rcpp::stop("");
  }
}

#endif
