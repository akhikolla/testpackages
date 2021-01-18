#' @useDynLib Immigrate
#' @importFrom Rcpp sourceCpp
NULL

#' IM4E
#'
#' This function performs IM4E(Iterative Margin-Maximization under Max-Min Entropy) algorithm.
#' @param xx model matrix of explanatory variables
#' @param yy label vector
#' @param epsilon criterion for stopping iteration, default to be 0.01
#' @param sig sigma used in algorithm, default to be 1
#' @param lambda lambda used in algorithm, default to be 1
#' @param max_iter maximum number of iteration 
#' @param removesmall whether remove features with small weights, default to be FALSE
#' @keywords IM4E
#' @return \item{w}{weight vector obtained by IM4E algorithm}
#' @return \item{iter_num}{number of iteration for convergence}
#' @return \item{final_c}{final cost value. Refer to the cost function in reference below for more details}
#' @export
#' @examples
#' data(park)
#' xx<-park$xx
#' yy<-park$yy
#' re<-IM4E(xx,yy)
#' print(re)
#' @references 
#' Bei Y, Hong P. Maximizing margin quality and quantity[C]//Machine Learning for Signal Processing (MLSP), 2015 IEEE 25th International Workshop on. IEEE, 2015: 1-6.


IM4E <- function(xx,yy,epsilon=0.01,sig=1,lambda=1,max_iter=10,removesmall=FALSE) {
  suppressWarnings(
  res<-(IM4ECpp(oneIM4E = one.IM4E,xx,yy,epsilon=0.01,
                 sig=1, lambda=1,max_iter=10,removesmall=FALSE)))
  class(res)<-"IM4E"
  return(res)
}
