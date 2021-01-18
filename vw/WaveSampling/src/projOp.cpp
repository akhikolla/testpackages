#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @encoding UTF-8
//' @title Projection operator
//'
//'
//' @description
//'
//' This operator projects the vector v orthogonally onto the line spanned by vector u.
//'
//' @param v vector projected.
//' @param u vector that define the line on which we project.
//' 
//' @details
//' 
//' The projection operator is defined by :
//' 
//' \deqn{proj_u(v) = \frac{\langle u , v \rangle}{\langle u, u \rangle} u}
//'  where \eqn{\langle . , . \rangle} is the inner product also written \eqn{u^\top v}.
//' 
//' @return The projection of the vector v onto the line spanned by the vector u.
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @references 
//' \url{https://en.wikipedia.org/wiki/Projection_(linear_algebra)}s
//' 
arma::vec projOp(arma::vec v,arma::vec u) {
  double num = arma::as_scalar(u.t()*v);
  double den = arma::as_scalar(u.t()*u);
  arma::vec scalar(u.size());
  scalar.fill(num/den);
  return scalar%u;
}


/*** R
u = c(0,1)
v = c(1,1)
projOp(v,u)
v - projOp(v,u)
*/
