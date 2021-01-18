#' Extreme Value random process generation. 
#' 
#' Extreme Value random process generation. 
#'
#' Generation of samples from unit Frechet processes.
#' @aliases rFrechet 
#' @param n Number of samples to generate.
#' @usage rFrechet(n) 
#' @examples
#' rFrechet(1000)
#' @export
 
rFrechet <- function(n){
  -1/log(runif(n))
}

#' Extreme Value random process generation. 
#' 
#' Extreme Value random process generation. 
#'
#' Generation of samples from  Max AR(theta) processes.
#' @aliases  rMaxAR
#' @param n Number of samples to generate.
#' @param theta Parameter of the MAX AR process, takes values between 0 and 1.
#' @usage rMaxAR(n,theta)
#' @examples
#' rMaxAR(1000,0.2)

#' @export

rMaxAR <- function(n,theta)
{
  x <- rFrechet(n)
  
  y <- rbind(rep(x[1],2),cbind((1-theta)*x[-n], x[-1]))
  
  y <- apply(y,1, max)
  
  y

}
