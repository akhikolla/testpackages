#include <RcppArmadillo.h>

#include "distUnitk.h"
#include "wpik.h"
#include "rowcolSumsiter.h"
#include "projOp.h"



// [[Rcpp::depends(RcppArmadillo)]]
//' @encoding UTF-8
//' @title Weakly associated vectors sampling
//'
//' @description
//'
//' Select a spread sample from inclusion probabilities using the weakly associated vectors sampling method.  
//'
//' @param X matrix representing the spatial coordinates. 
//' @param pik vector of the inclusion probabilities. The length should be equal to N.
//' @param bound a scalar representing the bound to reach. See Details. Default is 1.
//' @param tore an optional logical value, if we are considering the distance on a tore. See Details. Default is \code{TRUE}.
//' @param shift an optional logical value, if you would use a shift perturbation. See Details. Default is \code{FALSE}.
//' @param toreBound a numeric value that specify the size of the grid. Default is -1.
//' @param comment an optional logical value, indicating some informations during the execution. Default is \code{FALSE}.
//' @param fixedSize an optional logical value, if you would impose a fixed sample size. Default is \code{TRUE}
//'
//' @details
//' 
//' The main idea is derived from the cube method (Devill and Tillé, 2004). At each step, the inclusion probabilities vector \code{pik}
//' is randomly modified. This modification is carried out in a direction that best preserves the spreading of the sample.
//' 
//' A stratification matrix \eqn{\bf A} is computed from the spatial weights matrix calculated from the function \code{\link{wpik}}.
//' Depending if \eqn{\bf A} is full rank or not, the vector giving the direction is not selected in the same way.
//' 
//' If matrix \eqn{\bf A} is not full rank, a vector that is contained in the right null space is selected:
//' \deqn{ Null(\bf A) = \{ \bf x \in R^N | \bf A\bf x = \bf 0  \}. }
//' 
//' If matrix \eqn{\bf A} is full rank, we find \eqn{\bf v}, \eqn{\bf u} the singular vectors associated to the 
//' smallest singular value \eqn{\sigma } of \eqn{\bf A} such that
//' 
//' \deqn{ \bf A\bf v = \sigma \bf u,~~~~ \bf A^\top \bf u = \sigma \bf v.}
//' 
//' Vector \eqn{ \bf v } is then centered to ensure fixed sample size. At each step, inclusion probabilities is modified and at least on component is set to 0 or 1. Matrix \eqn{\bf A } is updated 
//' from the new inclusion probabilities. The whole procedure it repeated until it remains only one component that is not equal to 0 or 1.
//' 
//' For more informations on the options \code{tore} and \code{toreBound}, see \code{\link{distUnitk}}. If \code{tore} is set up \code{TRUE} and \code{toreBound} not specified the \code{toreBound} is equal to 
//' \deqn{N^{1/p}}
//' where \eqn{p} is equal to the number of column of the matrix \code{X}.
//' 
//' For more informations on the option \code{shift}, see \code{\link{wpik}}.
//' 
//' If \code{fixedSize} is equal \code{TRUE}, the weakest associated vector is centered at each step of the algorithm. This ensures that the size of the selected sample is equal to the sum of the inclusion probabilities.
//' 
//' @return A vector of size \eqn{N} with elements equal 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for non-chosen unit.
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @references 
//' Deville, J. C. and Tillé, Y. (2004). Efficient balanced sampling: the cube method. Biometrika, 91(4), 893-912
//' 
//' 
//' @seealso
//' \code{\link{wpik}}, \code{\link{distUnitk}}.
//' 
//' 
//' @examples
//' 
//' #------------
//' # Example 2D
//' #------------
//' 
//' N <- 50
//' n <- 15
//' X <- as.matrix(cbind(runif(N),runif(N)))
//' pik <- sampling::inclusionprobabilities(runif(N),n)
//' s <- wave(X,pik)
//' 
//' #------------
//' # Example 2D grid 
//' #------------
//' 
//' N <- 36 # 6 x 6 grid
//' n <- 12 # number of unit selected
//' x <- seq(1,sqrt(N),1)
//' X <- as.matrix(cbind(rep(x,times = sqrt(N)),rep(x,each = sqrt(N))))
//' pik <- rep(n/N,N)
//' s <- wave(X,pik, tore = TRUE,shift = FALSE)
//' 
//' #------------
//' # Example 1D
//' #------------
//' 
//' N <- 100
//' n <- 10
//' X <- as.matrix(seq(1,N,1))
//' pik <- rep(n/N,N)
//' s <- wave(X,pik,tore = TRUE,shift =FALSE,comment = TRUE)
//' 
//' 
//' @export
// [[Rcpp::export]]
arma::vec wave(const arma::mat& X,
               const arma::vec& pik,
               double bound = 1.0,
               bool tore = false,
               bool shift = false,
               double toreBound = -1,
               bool comment = false,
               bool fixedSize = true) {
  // INITIALIZE CONSTANT
  double la1 = 1e+200;
  double la2 = 1e+200;
  double la = 1e-9;
  double eps = 1e-13;
  int N = X.n_rows;
  
  
  //
  if(N >= 500){
    Rcpp::Rcout << "WARNING : Your population size is greater than 500 this could be quite time-consuming..." << std::endl;
  }
  
  //double is important in order to take 1/p
  double p = X.n_cols;
  // if the bound is specified.
  double tb = 0.0;
  if(toreBound == -1){
    tb = pow(N,1/p);
  }else{
    tb = toreBound;
  }

  // vector of indices of the number of dimension for extract matrix 
  arma::uvec ncol(p);
  for(int l = 0; l < p; l++){
    ncol(l) = l;
  }
  
  // vector of one
  arma::mat one = arma::ones<arma::mat>(N,1);
  
  // vector of pik that are gonna be updated
  arma::vec re(pik);
  arma::uvec i = arma::find(re > eps && re < (1-eps));
  
  unsigned int i_size = i.size();
  
  //INITIALIZING VARIABLE 
  arma::mat U;
  arma::vec s;
  arma::mat V;
  
  if(comment  == true){
    Rcpp::Rcout << "--- Sample selection ---" << std::endl;
  }
  
  // MAIN LOOP
  while(i_size > 1){
    
    if(comment  == true){
      if(i_size % 1 == 0){
        Rcpp::Rcout << "Current dimension of the Matrix W : " << i_size << std::endl;
      }
    }
    
    arma::mat A(i_size,i_size);
    one = arma::ones<arma::mat>(i_size,1);
    A = wpik(X.submat(i,ncol),re.elem(i),bound,tore,shift,tb);
    arma::mat D = diagmat(1/re.elem(i));
    A = A*D;
    
    
    // init vector for updating pik
    arma::vec u(i_size);
    
    /* QR OR SVD
     * 
     * svd_econ faster than a QR decomposition to extract null space.
     * 
     * TODO:
     * find a faster way to find the smallest associated
     * vector CHECK for LOPBCG and/or other methods.
     * 
     */
    arma::svd_econ(U,s,V,A,"right","dc");
    arma::uvec r = find(s >= eps);
    unsigned int rang = r.size();
    if(rang < i_size){
        // same for now but could be different in future dev.
        u = V.col(V.n_cols - 1);
    }else{
        u = V.col(V.n_cols - 1);
    }
    if(fixedSize == true){
      u = u - projOp(u,one);  
    }
    
    la1 = 1e+200;
    la2 = 1e+200;
    la = 1e-9;
    for(unsigned int k = 0; k < i.size(); k++){
      if(u[k]> 0){
        la1 = std::min(la1,(1.0-re[i[k]])/u[k]);
        la2 = std::min(la2,re[i[k]]/u[k]);
      }
      if(u[k]< 0){
        la1 = std::min(la1,-re[i[k]]/u[k]);
        la2 = std::min(la2,(re[i[k]]-1.0)/u[k]);
      }
    }
    if(Rcpp::runif(1)[0]<la2/(la1+la2)){
      la = la1;
    }else{
      la = -la2;
    }
    for(unsigned int k = 0; k < i.size(); k++){
      re[i[k]] = re[i[k]] + la*u[k];
      if(re[i[k]] < eps){
        re[i[k]] = 0;
      }
      if(re[i[k]] > (1-eps)){
        re[i[k]] = 1;
      }
    }
    i = arma::find(re > eps && re < (1-eps));
    i_size = i.size();
    
    //Check if the user stop the function (approximately every 10 iterations)
    if (i_size % 10 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    if(arma::sum(re.elem(i)) < (1-eps) && fixedSize == false){
      Rcpp::Rcout << "The algorithm end because the remaining sum of inclusion probabilites is equal to : " << arma::sum(re.elem(i)) << 
      "\nThe remaining inclusion probabilities are \n" << re.elem(i) <<
      "\nBernoulli distribution are used to select the units.\n";
      for(unsigned int tt = 0; tt < i.size(); tt++){
        re[i[tt]] = R::rbinom(1,re[i[tt]]);
      }
      Rcpp::Rcout << re.elem(i) << "\n";
      break;
    }

    
  }
  
  if(comment  == true){
    Rcpp::Rcout << re << std::endl;
    Rcpp::Rcout << "--- Sample selection finished ---" << std::endl;
  }
  
  
  return(arma::round(re));
}


/*** R


N <- 36
n <- 12
x <- seq(1,sqrt(N),1)
X <- as.matrix(cbind(runif(N),runif(N)))
X <- as.matrix(expand.grid(x,x))
pik <- rep(n/N,N)
# pik <- sampling::inclusionprobabilities(runif(N),n)
system.time(s <- wave(X,pik,tore = TRUE,shift = TRUE,fixedSize = TRUE))
sum(s)
plot(X)
points(X[s == 1,],pch = 16)



N <- 50
n <- 15
x <- as.matrix(runif(N),runif(N))
pik <- sampling::inclusionprobabilities(runif(N),n)
s <- wave(x,pik,comment = TRUE,fixedSize = FALSE)
sum(s)


############# EXAMPLE 1

N <- 36
n <- 12
x <- seq(1,sqrt(N),1)
X <- as.matrix(cbind(rep(x,times = sqrt(N)),rep(x,each = sqrt(N))))
pik <- rep(n/N,N)
s <- wave(X,pik,tore = T,shift =F,toreBound = -1, comment = TRUE)
s <- wave(X,pik,tore = T,shift =F,comment = TRUE)
plot(X)
points(X[s == 1,],pch = 16)

X <- as.matrix(cbind(runif(N),runif(N)))
s <- wave(X,pik,tore = F,shift =F,comment = TRUE,fixedSize = FALSE)
sum(s)
plot(X)
points(X[s == 1,],pch = 16)


N <- 225
n <- 75
x <- seq(1,sqrt(N),1)
X <- as.matrix(cbind(rep(x,times = sqrt(N)),rep(x,each = sqrt(N))))
pik <- rep(n/N,N)

# W <- wpik(X,pik,bound = 1,tore = TRUE,shift = TRUE,toreBound = 15 )
s <- wave(X,pik,tore = T,shift =T,comment = TRUE,fixedSize = FALSE)
plot(X)
points(X[s == 1,],pch = 16)


N <- 1000
n <- 100
# x <- seq(1,sqrt(N),1)
# X <- as.matrix(cbind(rep(x,times = sqrt(N)),rep(x,each = sqrt(N))))
X <- as.matrix(cbind(runif(N),runif(N)))
pik <- rep(n/N,N)
# W <- wpik(X,pik,bound = 1,tore = TRUE,shift = TRUE,toreBound = 15 )
system.time(s <- wave(X,pik,tore = T,shift =T,comment = T))
# WITH R open
# utilisateur     système      écoulé 
# 938.26        7.21      309.32 
# With R 3.6.1
# utilisateur     système      écoulé 
# 1073.12        6.43     1085.92
plot(X)
points(X[s == 1,],pch = 16)

  
N <- 25
n <- 5
x <- seq(1,sqrt(N),1)
X <- as.matrix(cbind(rep(x,times = sqrt(N)),rep(x,each = sqrt(N))))
pik <- rep(n/N,N)
s <- wave(X,pik,tore = TRUE,comment = TRUE)
s <- wave2(X,pik,tore = TRUE,comment = TRUE)
plot(X)
points(X[s == 1,],pch = 16)


rm(list = ls())
N <- 30
n <- 100
x <- seq(1,N,1)
y <- seq(1,N,1)
X <- as.matrix(expand.grid(x,y))

pik <- rep(n/(N*N),N*N)
W <- t(wpik(as.matrix(X),pik,bound = 1.0,tore = TRUE,shift = T,toreBound = N))
image(W)
system.time(test <- wave(as.matrix(X),pik, tore = TRUE,shift = T,comment = TRUE))
plot(X)
points(X[test ==1,],pch = 16)


#############
# 1D


N <- 100
n <- 10
x <- as.matrix(seq(1,N,1))
# x <- as.matrix(runif(N))
pik <- rep(n/N,N)
# pik <- sampling::inclusionprobabilities(runif(N),n)
s <- wave(x,pik,tore = T,shift =F,comment = T)
plot(x,rep(0,N))
points(x[s == 1,],rep(0,n),pch = 16)


*/
