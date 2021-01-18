#' Generate a row-wise stacked boolean partition matrix of JI rows and K columns
#'
#' @details The function generates a (JI x K) boolean partition matrix for I individuals, J stimuli and K categories. Note that J and K must be chosen so that J%%K = 0 (J must be a multiple of K).
#' This is an internal function, generally not to be called by the user.
#' @export 
#' @param I (integer) number of individuals
#' @param J (integer) number of trials
#' @param K (integer) number of levels for a categorical variables
#' @param type (character) method to generate the matrix: symmetric (default) or random
#' @return a (JI x K) boolean matrix
#' @examples

#' Z <- generate_Z(I = 2,J = 12,K = 4,type="symmetric")
#' print(Z)
#' 

generate_Z <- function(I,J,K,type=c("symmetric","random")){

  if(I<1 | J<1 | K<1)
    stop("Positive integers should be provided for I,J,K")
  if(K>J)
    stop("J must be equal or greater than K")
  if(J%%K>0)
    stop("The integer J should be a multiple of K (J%%K = 0)")
  type=match.arg(type)

  if(type=="random"){
    M <- matrix(stats::rnorm(I*J*K,sd=100),ncol=K)
    Z <- round(exp(M)/rowSums(exp(M)))
  }else{
    Z <- kronecker(kronecker(matrix(1,I,1),diag(K)),matrix(1,J/(K)))
  }

  return(Z)
}
