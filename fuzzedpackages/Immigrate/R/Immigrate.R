#' Immigrate
#'
#' This function performs IMMIGRATE(Iterative Max-Min Entropy Margin-Maximization with Interaction Terms ) algorithm.
#' IMMIGRATE is a hypothesis-margin based feature selection method with interaction terms.
#' Its weight matrix reflects the relative importance of features and their iteractions, which can be used for feature selection.  
#' @param xx model matrix of explanatory variables
#' @param yy label vector
#' @param w0 initial weight matrix, default to be diagonal matrix when missing
#' @param epsilon criterion for stopping iteration
#' @param sig sigma used in algorithm, default to be 1. Refer to the cost function in the link below for more details
#' @param max_iter maximum number of iteration 
#' @param removesmall whether to remove features with small weights, default to be FALSE
#' @param randomw0 whether to use randomly initial weights, default to be FALSE
#' @keywords Immigrate
#' @return \item{w}{weight matrix obtained by IMMIGRATE algorithm}
#' @return \item{iter_num}{number of iteration for convergence}
#' @return \item{final_c}{final cost value. Refer to the cost function in link below for more details}
#' @import Rcpp
#' @importFrom stats runif

#' @export
#' @examples
#' data(park)
#' xx<-park$xx
#' yy<-park$yy
#' re<-Immigrate(xx,yy)
#' print(re)
#' @references Zhao, Ruzhang, Pengyu Hong, and Jun S. Liu. "IMMIGRATE: A Margin-based Feature Selection Method with Interaction Terms." Entropy 22.3 (2020): 291.
#' @seealso Please refer to \url{https://www.mdpi.com/1099-4300/22/3/291/htm} for more details.
#' @seealso Please refer to \url{https://github.com/RuzhangZhao/Immigrate/} for implementation demo.
Immigrate<-function(xx,yy,w0,epsilon=0.01,sig=1, max_iter=10,removesmall=FALSE,randomw0=FALSE){
  p<-ncol(xx)
  if (missing(w0)){
    if (randomw0){
      w0 <-matrix(runif(p*p),p,p)
    }else{
      w0 <-diag(1,p,p)
    }
    w0<-w0/sqrt(sum(w0^2))
  }
  suppressWarnings(
  res<-ImmigrateCpp(oneImmigrate=one.Immigrate,xx,yy,w0 = w0,epsilon=epsilon,
                      sig=sig, max_iter=max_iter,removesmall=removesmall))
  class(res)<-"Immigrate"
  return(res)
}
