#include <RcppArmadillo.h>
#include "distUnitk.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @encoding UTF-8
//' @title Stratification matrix from inverse inclusion probabilities
//'
//' @description
//'
//' The stratification matrix is calculated from the inverse inclusion probabilities. It is a direct
//' implementation of the spatial weights specified in Tillé et al., (2018).
//'
//' @param X matrix representing the spatial coordinates. 
//' @param pik vector of the inclusion probabilities. The length should be equal to N.
//' @param tore an optional logical value, if we are considering the distance on a tore. Default is \code{FALSE}.
//' @param shift an optional logical value, if you would use a shift perturbation. See Details for more informations. Default is \code{FALSE}.
//' @param toreBound a numeric value that specify the size of the grid. Default is -1.
//' 
//' @details
//' 
//' Entries of the stratification matrix indicates how the units are close from each others. Hence a large value \eqn{w_{kl}} means that the unit \eqn{k} 
//' is close to the unit \eqn{l}. This function considers that if unit \eqn{k} were selected in the sample drawn from the population then
//' \eqn{k} would represent \eqn{1/\pi_k} units in the population and, as a consequence, it would be natural to consider that
//' \eqn{k} has \eqn{n_k =  (1/\pi_k - 1)} neighbours in the population. The \eqn{n_k} neighbours are the nearest neighbours of \eqn{k} according to distances.
//' The weights are so calculated as follows :
//' 
//' \itemize{
//'   \item \eqn{ w_{kl} = 1} if unit \eqn{l \in N_{\lfloor n_k \rfloor }}
//'   \item \eqn{ w_{kl} = n_k - \lfloor n_k \rfloor} if unit \eqn{l} is the \eqn{\lceil n_k \rceil} nearest neighbour of \eqn{k}.
//'   \item \eqn{w_{kl} = 0} otherwise.
//' }
//' 
//' 
//' \eqn{\lfloor n_k \rfloor } and \eqn{\lceil n_k \rceil} are the inferior and the superior integers of \eqn{n_k}.
//' 
//' The option \code{shift} add a small normally distributed perturbation \code{rnorm(0,0.01)} to the coordinates
//' of the centroid of the stratum considered. This could be useful if there are many unit that have the same distances.
//' Indeed, if two units have the same distance and are the last unit before that the bound is reached, then the weights
//' of the both units is updated. If a shift perturbation is used then all the distances are differents and only one unit
//' weight is update such that the bound is reached. 
//' 
//' The shift perturbation is generated at the beginning of the procedure such that each stratum is shifted by the same perturbation.
//'
//' 
//' @return A sparse matrix representing the spatial weights.
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @references 
//' Tillé, Y., Dickson, M.M., Espa, G., and Guiliani, D. (2018). Measuring the spatial balance of a sample: A new measure based on Moran's I index.
//' \emph{Spatial Statistics}, 23, 182-192. \url{https://doi.org/10.1016/j.spasta.2018.02.001}
//' 
//' @seealso
//' \code{\link{wpik}}, \code{\link{distUnitk}}, \code{\link{wave}}.
//' @examples
//' 
//' N <- 25
//' n <- 5
//' X <- as.matrix(cbind(runif(N),runif(N)))
//' pik <- sampling::inclusionprobabilities(runif(N),n)
//' W <- wpikInv(X,pik)
//' 
//' @export
// [[Rcpp::export]]
arma::sp_mat wpikInv(arma::mat X,
                  arma::vec pik,
                  bool tore = false,
                  bool shift = false,
                  double toreBound = -1){

  /*
  * Initializing variable
  */
  int N = X.n_rows;
  arma::sp_mat W = arma::sp_mat(N,N);

  /*
  * If jiiter is is equal true then we
  *  add a perturbation on the central element of each strata
  *  specified here because we need to ensure that we add the
  *  same perturbation to each points.
  */
  double tb1(0.0);
  double tb2(0.0);
  if(shift == true){
    tb1 = R::rnorm(0.0, 0.01);
    tb2 = R::rnorm(0.0, 0.01);
  }
  /*
  * Main loop on the row of the Matrix X
  */
  for(int k = 1; k <= N; k++){
    arma::vec d(N);

    /*
    * Choose between the distance
    * -----------------------a tore.
    * --------------- add a jiter.
    */

    if(shift == true){
      double tmp1 = X(k-1,0);
      double tmp2 = X(k-1,1);
      X(k-1,0) = X(k-1,0) + tb1;
      X(k-1,1) = X(k-1,1) + tb2;
      if(tore == true){
        d = distUnitk(X,k,true,toreBound);
      }else{
        d = distUnitk(X,k,false,toreBound);
      }

      X(k-1,0) = tmp1;
      X(k-1,1) = tmp2;
    }else{
      if(tore == true){
        d = distUnitk(X,k,true,toreBound);
      }else{
        d = distUnitk(X,k,false,toreBound);
      }
    }

    
    /*
    * idx is the stable sort index of the distance.
    * such that d(idx) = 0 1 1 4 4 5 5 ....
    * CAREFUL : index of the distance 0-based.
    */
    // arma::uvec idx =  arma::stable_sort_index(d); // sort the distance with exact index even if we have ties

    arma::uvec d_unique = arma::find_unique(d,true);
    arma::vec d_unique_sorted  = arma::sort(d.elem(d_unique));
    arma::vec w(N,arma::fill::zeros); // returned vec

    double k_i = 1/pik[k-1];
    int k_i_floor = floor(1/pik[k-1]);
    // int k_i_ceil = ceil(1/pik[k-1]);

    arma::uvec modif;
    arma::vec w_tmp(N);
    int tt = 0;
    int nbr_element = 0;
    do{
      w_tmp.fill(0.0);
      // modif = find(d == d[d_unique[tt]]);
      modif = arma::find(d == d_unique_sorted[tt]);
      w.elem(modif) += 1.0;
      w_tmp.elem(modif) +=1.0;
      tt = tt + 1;
      nbr_element += modif.size();
    } while (nbr_element <= k_i_floor && nbr_element < N);

    double modif_value = (k_i-sum(w-w_tmp))/sum(w_tmp);

    w.elem(modif) -= 1.0;
    w.elem(modif) += modif_value;
    w[k-1] = 0.0;

    /*
    * column major implementation seems to be faster to do like that instead of adding rows.
    */
    W.col(k-1) = w;
  }
  return W.t();
}



/*** R

X <- cbind(runif(1000),runif(1000))
pik <- sampling::inclusionprobabilities(runif(1000),100)
d <- array(rep(0,1000*1000),c(1000,1000))
for(i in 1:1000){
  d[i,] <- distUnitk(X,k =i,tore = FALSE,toreBound = 0)
}


rm(list = ls())

wpik2R <- function(d,pik)
{
  N=length(pik)
  w=array(0,c(N,N))
  for(i in 1:N)
  {
    rr=rank(d[i,],ties.method="min")
    ww=as.integer(rr<=min(ceiling(1/pik[i]),N))
    dec=as.integer(rr==max(rr*ww))
    if(sum(dec)>0)  ww[dec==1]=ww[dec==1]*(1/pik[i]-sum(ww-dec))/sum(dec)
    w[i,]=ww
    diag(w)=0
  }
  w
}

x <- seq(1,5,1)
y <- seq(1,5,1)
X <- as.matrix(expand.grid(x,y))
pik <- rep(0.1352345,25)
d <- array(rep(0,5*5),c(25,25))
for(i in 1:25){
  d[i,] <- distUnitk(X,k =i,tore = TRUE,toreBound = 5)
}

as(wpik2R(d,pik),"sparseMatrix") == wpikInv(X,pik = pik,tore = TRUE,shift = FALSE,toreBound = 5)


X <- cbind(runif(25),runif(25))
pik <- sampling::inclusionprobabilities(runif(25),5)
pik[1] <- 0.015345
d <- array(rep(0,5*5),c(25,25))
for(i in 1:25){
  d[i,] <- distUnitk(X,k =i,tore = FALSE,toreBound = 0)
}

as(wpik2R(d,pik), "sparseMatrix") == wpikInv(X,pik = pik,tore = FALSE,shift = FALSE,toreBound = 0)
image(wpikInv(X,pik = pik,tore = FALSE,shift = FALSE,toreBound = 0))


### time

X <- cbind(runif(1000),runif(1000))
pik <- sampling::inclusionprobabilities(runif(1000),100)
d <- array(rep(0,1000*1000),c(1000,1000))
for(i in 1:25){
  d[i,] <- distUnitk(X,k =i,tore = FALSE,toreBound = 0)
}


system.time(as(wpik2R(d,pik), "sparseMatrix"))
# utilisateur     système      écoulé 
# 12.00        2.18       14.18 
system.time(W <- wpikInv(X,pik = pik,tore = FALSE,shift = FALSE,toreBound =0))
# utilisateur     système      écoulé 
# 0.81        0.00        0.81 





*/

