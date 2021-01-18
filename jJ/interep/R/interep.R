#' @useDynLib interep, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' fit generalized estimaitng equations with given tuning parameters
#'
#' This function makes predictions for generalized estimating equation with a given value of lambda.
#' Typical usage is to have the cv.interep function compute the optimal lambda, then provide it to
#' the interep function.
#' @importFrom stats gaussian
#' @importFrom MASS ginv
#' @param e matrix of environment factors.
#' @param g matrix of omics factors. In the case study, the omics measurements are lipidomics data.
#' @param y the longitudinal response.
#' @param beta0 the inital coefficient vector.
#' @param corre the working correlation structure that is used in the estimation algorithm. interep provides three choices for the
#' working correlation structure: "a" as AR-1", "i" as "independence" and "e" as "exchangeable".
#' @param pmethod the penalization method. "mixed" refers to MCP penalty to individual main effects and group MCP penalty to interactions; "individual" means MCP penalty to all effects.
#' @param lam1 the tuning parameter lambda1 for individual predictors.
#' @param lam2 the tuning parameter lambda2 for interactions.
#' @param maxits the maximum number of iterations that is used in the estimation algorithm.  The default value is 30
#' @return
#' \item{coef}{the coefficient vector.}
#' @references
#' Zhou, F., Ren, J., Li, G., Jiang, Y., Li, X., Wang, W.and Wu, C. (2019). Penalized variable selection for Lipid--environment interactions in a longitudinal lipidomics study.
#' \href{https://www.mdpi.com/2073-4425/10/12/1002/htm}{\emph{Genes}, 10(12), 1002}
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S. and Wu, C. (2020) Gene–Environment Interaction: a Variable Selection Perspective.
#' \href{https://arxiv.org/abs/2003.02930}{\emph{Epistasis}, Methods in Molecular Biology. Humana Press. (Accepted)}
#'
#' @examples
#' data("dat")
#' e=dat$e
#' g=dat$z
#' y=dat$y
#' beta0=dat$coef
#' index=dat$index
#' b = interep(e, g, y,beta0,corre="e",pmethod="mixed",lam1=dat$lam1, lam2=dat$lam2,maxits=30)
#' b[abs(b)<0.05]=0
#' pos = which(b != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' @export

interep <- function(e,g,y,beta0,corre,pmethod,lam1,lam2,maxits){
  q=dim(e)[2]
  n=dim(y)[1]
  k=dim(y)[2]
  p1=dim(g)[2]

  x=cbind(e,g)
  for (i in 1:p1) {
    for (j in 1:q) {
      x=cbind(x,e[,j]*g[,i])
    }
  }

  x=scale(x)

  #==========================================reformat the data===============================#
  data=reformat(k,y,x)
  y=data$y
  x=data$x

  p=dim(x)[2]
  k=rep(k,n)

  converge=F
  iter=0
  beta.new=beta0
  #=========================================PQIF======================================#
  while ((!converge) & (iter < maxits)) {
    beta = beta.new
    Score=ScoreU(n,k,y,x,p,beta,corre)
    U=Score$U
    dU=Score$dU

    E=penalty(x,n,p,q,beta,lam1,pmethod,p1,lam2)

    mat=dU + n*E
    mat[abs(mat)<0.000001]=0
    beta.new = beta + MASS::ginv(mat)%*%(U - n*E%*%beta)

    diff=mean(abs(beta.new-beta))
    converge = (diff < 1e-3)
    iter = iter+1
    cat("iter",iter,"diff",diff,"\n")
  }
  coef=beta.new
  return(coef)
}
