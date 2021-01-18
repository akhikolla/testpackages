#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Squared Euclidean distances of the unit k.
//'
//' @description
//' Calculate the squared Euclidean distance from unit \eqn{k} to the other units.
//' 
//'
//' @param X matrix representing the spatial coordinates. 
//' @param k the unit index to be used.
//' @param tore an optional logical value, if we are considering the distance on a tore. See Details.
//' @param toreBound an optional numeric value that specify the length of the tore.
//'
//'
//' @details
//' 
//' Let \eqn{\mathbf{x}_k,\mathbf{x}_l} be the spatial coordinates of the unit \eqn{k,l \in U}. The classical euclidean distance is given by
//' 
//' \deqn{d^2(k,l) = (\mathbf{x}_k - \mathbf{x}_l)^\top (\mathbf{x}_k - \mathbf{x}_l). }
//' 
//' When the points are distributed on a \eqn{N_1 \times N_2} regular grid of \eqn{R^2}.
//' It is possible to consider the units like they were placed on a tore. It can be illustrated by Pac-Man passing through the wall to get away from ghosts. Specifically,
//' we could consider two units on the same column (resp. row) that are on the opposite have a small distance,
//' 
//' \deqn{ d^2_T(k,l) = min( (x_{k_1} - x_{l_1})^2,
//'                       (x_{k_1} + N_1 - x_{l_1})^2,
//'                       (x_{k_1} - N_1 - x_{l_1})^2) +}
//' \deqn{ min( (x_{k_2} - x_{l_2})^2,
//'                       (x_{k_2} + N_2 - x_{l_2})^2,
//'                       (x_{k_2} - N_2 - x_{l_2})^2).}
//'
//' The option \code{toreBound} specify the length of the tore in the case of \eqn{N_1 = N_2 = N}. 
//' It is omitted if the \code{tore} option is equal to \code{FALSE}.
//'
//' @return a vector of length \eqn{N} that contains the distances from the unit \eqn{k} to all other units.
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' 
//' @seealso
//' \code{\link{wpik}}, \code{\link{wave}} and \code{\link[stats]{dist}}.
//'
//' @examples
//'   N <- 5
//'   x <- seq(1,N,1)
//'   X <- as.matrix(expand.grid(x,x))
//'   distUnitk(X,k = 2,tore = TRUE,toreBound = 5)
//'   distUnitk(X,k = 2,tore = FALSE,toreBound = -1)
//' @export
// [[Rcpp::export]]
arma::vec distUnitk(arma::mat X,
                    int k,
                    bool tore,
                    double toreBound){
  
  //initializing variables
  unsigned int N = X.n_rows;
  unsigned int p = X.n_cols;
  
  arma::vec unitkX(N);
  arma::vec unitkY(N);
  
  
  
  // Matrix with coordinates equal to the unit k
  arma::vec tmp(N);
  arma::mat coord_unitk(N,p);
  for(unsigned int j = 0; j < p; j++){
    coord_unitk.col(j) = tmp.fill(X((k-1),j));
  }
  
  
  //out vector
  arma::vec dist(N);
  
  // temporary matrix for the tore option
  arma::mat x1(N,p);
  arma::mat x2(N,p);
  arma::mat x3(N,p);
  arma::mat x(N,p);
  if(tore == true){
    
    // toreBound vector
    arma::vec toreBound_vec(N);
    toreBound_vec.fill(toreBound);

    // calculate the three distance
    for(unsigned int k = 0; k < p; k++){
      x1.col(k) = (coord_unitk.col(k)  - X.col(k))% (coord_unitk.col(k)  - X.col(k));
      x2.col(k) = (coord_unitk.col(k) + toreBound_vec - X.col(k))%(coord_unitk.col(k) + toreBound_vec - X.col(k));
      x3.col(k) = (-coord_unitk.col(k) + toreBound_vec + X.col(k))%(-coord_unitk.col(k) + toreBound_vec + X.col(k));
    }
    
    // choose the minimal one
    for(unsigned int i = 0;i < N; i++){
      for(unsigned int j = 0;j < p; j++){
        x(i,j) = std::min(x1(i,j),x2(i,j));
        x(i,j) = std::min(x(i,j),x3(i,j));
      }
    }
    
    //sum over the possible dimension
    dist = arma::sum(x,1);

  }else{
  
    // Euclidean distance without tore so just on each dimension
    for(unsigned int k = 0; k < p; k++){
      x1.col(k) = (coord_unitk.col(k)  - X.col(k))% (coord_unitk.col(k)  - X.col(k));
    }
    
    // sum over the possible dimension
    dist = arma::sum(x1,1);
  }
  
  return(dist);
}




/*** R
x <- seq(1,2,1)
y <- seq(1,2,1)
X <- as.matrix(expand.grid(x,y))
distUnitk(X,1,tore = TRUE,toreBound = NA)



N <- 5
x <- seq(1,N,1)
X <- as.matrix(expand.grid(x,x))
distUnitk(X,2,tore = TRUE,toreBound = 5)
distUnitk(X,2,tore = TRUE,toreBound = 5)


X <- as.matrix(seq(1,10,1))
X <- as.matrix(runif(10))
dist(X)
distUnitk(X,1,tore = TRUE,toreBound = 10)
distUnitk2(X,4,tore = FALSE,toreBound = NA)




N <- 100
x <- seq(1,N,1)
X <- as.matrix(expand.grid(x,x))
system.time(test1 <- distUnitk(X,100,tore = TRUE,toreBound = 100))
system.time(test2 <- distUnitk2(X,2,tore = TRUE,toreBound = 100))


*/
