#include <RcppArmadillo.h>
#include "rowcolSumsiter.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @encoding UTF-8
//' @title Spreading measure based on Moran's \eqn{I} index
//'
//' @description
//'
//' This function implements the spreading measure based on Moran's \eqn{I} index.
//'
//' @param W a stratification matrix inheriting from \code{\link[Matrix]{sparseMatrix}} that represents the spatial weights. See \code{\link{wpik}}.
//' @param s a vector of size \eqn{N} with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for non-chosen units.
//'
//'
//' @details
//' 
//' This index is developped by Tillé et al. (2018) and measure the spreading of a sample drawn from a population.
//' It uses a corrected version of the traditional Moran's \eqn{I} index. Each row of the matrix \eqn{\bf W} should represents a stratum. Each 
//' stratum is defined by a particular unit and its neighbouring units. See \code{\link{wpik}}.
//' The spatial balance measure is equal to
//' 
//' \deqn{I_B =\frac{( \bf s- \bar{s}_w)^\top  W ( s- \bar{s}_w)}{\bf \sqrt{( s- \bar{s}_w)^\top  D ( s- \bar{s}_w) ( s- \bar{s}_w)^\top  B ( s- \bar{s}_w)}},}
//' 
//' where \eqn{\bf D} is the diagonal matrix containing the \eqn{w_i}, 
//' 
//' \deqn{ \bf \bar{s}_w =  1 \frac{ s^\top  W  1}{ 1^\top  W  1}}
//' 
//' and 
//' 
//' \deqn{ \bf B =  W^\top  D^{-1}  W - \frac{ W^\top  1 1^\top  W}{1^\top  W  1}.}
//' 
//' To specify the spatial weights uses the argument \code{W}.
//'
//' @return A numeric value that represents the spatial balance. It could be any real value between -1 (spread) and 1 (clustered). 
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @references
//' Tillé, Y., Dickson, M.M., Espa, G., and Guiliani, D. (2018). Measuring the spatial balance of a sample: A new measure based on Moran's I index.
//' \emph{Spatial Statistics}, 23, 182-192. \url{https://doi.org/10.1016/j.spasta.2018.02.001}
//' 
//' @seealso
//' \code{\link{wpik}}
//' 
//' @examples
//'   N <- 36
//'   n <- 12
//'   x <- seq(1,sqrt(N),1)
//'   X <- expand.grid(x,x)
//'   pik <- rep(n/N,N)
//'   W <- wpik(as.matrix(X),pik,bound = 1,tore = TRUE,shift = FALSE,toreBound = sqrt(N))
//'   W <- W - diag(diag(W))
//'   s <- wave(as.matrix(X),pik,tore = TRUE,shift = TRUE,comment = TRUE)
//'   IB(W,s)
//' 
//' @export
// [[Rcpp::export]]
double IB(const arma::sp_mat& W,const arma::vec& s){
  
  // initializing temp variables
  int N = s.size();
  
  // rowsum of the matrix W
  arma::vec wi = rowSumsiter(W);
  
  // weighted mean sample
  arma::vec sb(N);
  // total sum of the matrix 
  double tot =0.0;
  for(arma::sp_mat::const_iterator it = W.begin(); it != W.end(); ++it) {
    tot = tot + *it;
  }
  
  sb.fill(arma::sum(wi%s)/tot);
  
  // centered sample
  arma::vec z = s-sb;
  
  // unity vector
  arma::vec u(N);
  u.fill(1.0);
  
  //matrix B
  arma::mat B;
  B = W.t()*diagmat(1/wi)*W - (W.t()*u)*u.t()*W/tot;
  
  // temporary values
  double num = as_scalar(z.t()*W*z);
  double den1 = as_scalar(z.t()*diagmat(wi)*z);
  double den2 = as_scalar(z.t()*B*z);
  double den = sqrt(den1*den2);
  
  double I = num/den;
  return I;
}

/*** R

N <- 576
n <- 100
x <- seq(1,sqrt(N),1)
X <- expand.grid(x,x)
pik <- rep(n/N,N)
W <- wpik(as.matrix(X),pik,bound = 1,tore = TRUE,shift = FALSE,toreBound = sqrt(N))
W <- W - diag(diag(W))
# utilisateur     système      écoulé 
# 108.59        0.70      109.30
system.time(s <- wave(as.matrix(X),pik,tore = TRUE,shift = TRUE,comment = TRUE))

system.time(I1 <- IB(W,s))
plot(X)
points(X[s == 1,],pch = 16)

*/
