#' L_p q-Wasserstein Distance
#' 
#' Compute the q-Wasserstein distance between persistence diagrams using an arbitrary L_p norm 
#' as ground metric.
#'
#' @param d1 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param d2 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param dimension Dimension of the topological features of interest (0 for connected components, 1 for cycles etc).
#' @param q Order of the q-Wasserstein distance.
#' @param p Order of the L_p norm to be used as a ground metric in the computation of the Wasserstein distance.
#' @return The value for the L_p q-Wassesterstein between \code{d1} and \code{d2}.
#' @author Tullia Padellini, Francesco Palini. The included C++ library is authored by Michael Kerber, Dmitriy Morozov, and Arnur Nigmetov
#' @details This function provides an R interface for the efficient C++ library `HERA` by Michael Kerber, Dmitriy Morozov, and Arnur Nigmetov (\url{https://bitbucket.org/grey_narn/hera/src/master/}).
#' @examples 
#' diag1 <- matrix(c(1,1,1,0,2,3,2,2.5,4), ncol = 3, byrow = FALSE)
#' diag2 <- matrix(c(1,1,0,1,1,2), ncol = 3, byrow = FALSE)
#' wasserstein.distance(diag1, diag2, dimension = 1, q = 1, p = 2)
#' @references
#' \insertRef{kerber2017geometry}{kernelTDA}
#' @export
wasserstein.distance = function(d1, d2, dimension,  q, p = 2) {
  
  
  d1 = d1[d1[,1]==dimension, 2:3, drop = F]
  d2 = d2[d2[,1]==dimension, 2:3, drop = F]
  
  
  wasserstein_distance(Diag1 = d1, Diag2 = d2, q = q, internal_p = p, 
                       delta = 0.01, initial_eps = 1.0, eps_factor = 5.0)
  
}
