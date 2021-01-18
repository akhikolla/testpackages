#ifndef RIEMFUNC_INVEQUIV_H
#define RIEMFUNC_INVEQUIV_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "../include/riemfactory.hpp"

///////////////////////////////////////////////////////////////
inline arma::mat riemfunc_invequiv(arma::vec x, int m, int n, std::string name){
  if (name=="euclidean"){
    return(euclidean_invequiv(x,m,n));
  } else if (name=="sphere"){
    return(sphere_invequiv(x,m,n));
  } else if (name=="spd"){
    return(spd_invequiv(x,m,n));
  } else if (name=="grassmann"){
    return(grassmann_invequiv(x,m,n));
  } else if (name=="stiefel"){
    return(stiefel_invequiv(x,m,n));
  } else {
    Rcpp::Rcout << "RiemBase::riemfunc_invequiv : " << name << " is not yet implemented." << std::endl;
    Rcpp::stop("");  
  }
}
#endif
