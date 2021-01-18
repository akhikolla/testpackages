#' Compute the matrix of distances D for kappa parameters
#'
#' @details The function compute the distance of the Y-trajectories from the distractor and target points. This is an internal function, generally not to be called by the user.
#' @export
#' @param Y (matrix) N x JI matrix of observed trajectories
#' @param y_T (numeric) position in angles of the target
#' @param y_D (numeric) position in angles of the distractor
#' @return a N x JI matrix containing the delta values for each data point in Y
#' @examples 
#' 
#' ## Generate a generic matrix Y of I = 5 individuals and J = 1 trajectories (N = 61)
#' I <- 5; N <- 61
#' y_T <- pi/4; y_D <- (3*pi)/4
#' Y <- matrix(stats::rnorm(n = N*I,mean = (y_T+y_D)/2,sd = 10),N,I)
#' DY <- compute_D(Y=Y,y_T=y_T,y_D=y_D)
#' 


compute_D <- function(Y=NULL,y_T=pi/4,y_D=(3*pi)/4){
  
  DY <- mapply(function(i){ifelse(Y[,i] >= ((y_T+y_D)/2),
                            abs(Y[,i]-y_D),
                            abs(Y[,i]-y_T))},1:dim(Y)[2])
  
  return(DY)
}
