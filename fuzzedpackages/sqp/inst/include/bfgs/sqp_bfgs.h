#ifndef __sqp_bfgs_bfgs_included__        // if include guard for '/bfgs/sqp_bfgs.h' is undefined
#define __sqp_bfgs_bfgs_included__        // define include guard for '/bfgs/sqp_bfgs.h'

#include "sqp.h"

namespace sqp {
namespace bfgs{



inline void bfgs_update(arma::mat &hessian,
                        arma::vec &old_y,
                        arma::vec &new_y,
                        arma::vec &old_gradient,
                        arma::vec &new_gradient,
                        const bool constraint_adjustment = true)
{
  
  const arma::vec s = (new_y) - (old_y);
  const arma::vec y = (new_gradient) - (old_gradient);
  
  const double sy = arma::as_scalar(s.t()*y);
  
  const arma::vec Hs = (hessian)*s;
  
  const double sHs = arma::as_scalar(s.t()*(hessian)*s);
  
  const arma::vec* eta;
  
  if(sy < 0.2*sHs && constraint_adjustment)
  {
    const double theta = (0.8*sHs)/(sHs - sy);
    const arma::vec eta_val = theta * y + (1-theta)*Hs;
    eta = &eta_val;
    
  } else
  {
    eta =  &y;
  }
  
  (hessian) += ((*eta)*(*eta).t())/arma::as_scalar(s.t() * (*eta) ) - (Hs*s.t()*(hessian))/sHs;
  
}



}
}

#endif                                   // end of include guard for '/bfgs/sqp_bfgs.h'


