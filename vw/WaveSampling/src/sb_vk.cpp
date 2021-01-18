#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @encoding UTF-8
//' @title Values \eqn{v_k} to compute the Spatial balance
//' 
//' @description
//' 
//' Calculates the \eqn{v_k} values of the spatial balance developped by Stevens and Olsen (2004) and suggested by Grafström et al. (2012).
//' 
//' @param pik vector of the inclusion probabilities. The length should be equal to \eqn{N}.
//' @param X matrix representing the spatial coordinates.
//' @param s A vector of size \eqn{N} with elements equal 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for non-chosen unit.
//' 
//' @details
//' 
//' The spatial balance measure based on the Voronoï polygons is defined by 
//' 
//' \deqn{B(S) = \frac{1}{n}\sum_{k\in S} (v_k -1)^2 .}
//' 
//' The function return the \eqn{v_k} values and is mainly based on the function \code{\link[BalancedSampling:sb]{sb}} of the package \code{BalancedSampling} Grafström and Lisic (2019).
//' 
//' @return A vector of size \eqn{N} with elements equal to the \eqn{v_k} values. If the unit is not selected then the value is equal to 0.
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @references 
//' 
//' Grafström, A., Lundström, N.L.P. and Schelin, L. (2012). Spatially balanced sampling through the Pivotal method. 
//' \emph{Biometrics}, 68(2), 514-520
//' 
//' Grafström, A., Lisic J. (2019). BalancedSampling: Balanced and Spatially Balanced Sampling. R package version 1.5.5.
//' https://CRAN.R-project.org/package=BalancedSampling
//' 
//' Stevens, D. L. Jr. and Olsen, A. R. (2004). Spatially balanced sampling of natural resources.
//' \emph{Journal of the American Statistical Association 99, 262-278}
//' 
//' @seealso
//' \code{\link[BalancedSampling:sb]{BalancedSampling::sb}}
//' 
//' 
//' @examples
//' N <- 50
//' n <- 10
//' X <- as.matrix(cbind(runif(N),runif(N)))
//' pik <- sampling::inclusionprobabilities(runif(N),n)
//' s <- wave(X,pik)
//' v <- sb_vk(pik,X,s)
//' 1/n*sum((v[which(v != 0)]-1)^2)
//' BalancedSampling::sb(pik,X,which(s == 1))
//' 
//' @export
// [[Rcpp::export]]
arma::vec sb_vk(arma::vec pik,arma::mat X,arma::vec s) {
  
  int N = X.n_rows;
  int ncol = X.n_cols;
  
  
  arma::uvec s1 = find(s == 1); // 0 based
  int n = s1.n_elem;
  
  arma::vec vk(N,arma::fill::zeros);
  arma::vec nearest(n,arma::fill::zeros);
  int nrNearest = 0;
  double mindist = 1e+200;
  double d = 0.0;
  
  // loop on the population
  for(int i=0;i<N;i++){ 
    nrNearest = 0; 
    mindist = 1e+200;
    for(int j=0;j<n;j++){ 
      d = 0.0;
      // loop on the number of coordinates 
      for(int k=0;k<ncol;k++){ 
        // distance i to the jth element
        d += pow(X(i,k)-X(s1(j),k),2); 
      }
      if(d==mindist){
        nearest(nrNearest) =s1(j);
        nrNearest++;
      }
      if(d<mindist){
        mindist = d;
        nearest(0) = s1(j);
        nrNearest = 1;
      }
    }
    // loop on the number of point in the Voronoi polygons
    for(int k=0;k<nrNearest;k++){ 
      vk(nearest(k)) = vk(nearest(k)) + pik(i)/nrNearest;
    }
  }
  
  return(vk);
  
}


/*** R

X <- as.matrix(cbind(runif(50),runif(50)))
pik <- rep(10/50,50)
pik <- inclusionprobabilities(runif(50),10)
s <- wave(X,pik)

plot(X)
points(X[s == 1,],pch = 16)

v <- sb_vk(pik,X,s)
1/10*sum((v[which(v != 0)]-1)^2)

BalancedSampling::sb(pik,X,which(s == 1))

*/

