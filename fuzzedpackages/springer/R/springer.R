#' @useDynLib springer, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' fit the model with given tuning parameters
#'
#' This function performs penalized variable selection for longitudinal data based on generalized estimating equation (GEE) or quadratic inference functions (QIF) with a given value of lambda.
#' Typical usage is to first obtain the optimal lambda using cross validation, then provide it to the springer function.
#' @importFrom stats gaussian
#' @importFrom MASS ginv
#' @param clin a matrix of clinical covariates. The default value is NULL. Whether to include the clinical covariates is decided by user.
#' @param e a matrix of environment factors.
#' @param g a matrix of genetic factors.
#' @param y the longitudinal response.
#' @param beta0 the initial coefficient vector
#' @param func the framework to obtain the score equation.  Two choices are available: "GEE" and "QIF".
#' @param corr the working correlation structure adopted in the estimation algorithm. The springer provides three choices for the
#' working correlation structure: exchangeable, AR-1,and independence.
#' @param structure Three choices are available for structured variable selection. "bilevel" for sparse-group selection on both group-level and individual-level. "group" for selection on group-level only. "individual" for selection on individual-level only.
#' @param lam1 the tuning parameter \eqn{\lambda_1} for individual-level penalty applied to genetic factors.
#' @param lam2 the tuning parameter \eqn{\lambda_2} for group-level penalty applied to gene-environment interactions.
#' @param maxits the maximum number of iterations that is used in the estimation algorithm. The default value is 30.
#' @param tol The tolerance level. Coefficients with absolute values that are smaller than the tolerance level will be set to zero. The adhoc value can be chosen as 0.001.
#' @return
#' \item{coef}{the coefficient vector.}
#'
#' @details Look back to the data model described in "\code{\link{dat}}":
#'  \deqn{Y_{ij}= \alpha_0 + \sum_{m=1}^{t}\theta_m Clin_{ijm} + \sum_{u=1}^{q}\alpha_u E_{iju} + \sum_{v=1}^{p}\eta_v^\top Z_{ijv}+\epsilon_{ij},}
#' where \eqn{Z_{ijv}} contains the \eqn{v}th genetic main factor and its interactions with the \eqn{q} environment factors for the \eqn{j}th measurement on the \eqn{i}th subject
#' and \eqn{\eta_{v}} is the corresponding coefficient vector of length \eqn{1+q}.
#'
#' When {structure="bilevel"}, variable selection for genetic main effects and gene-environment interactions under the longitudinal response will be conducted on both individual and group levels (bi-level selection):
#' \itemize{
#' \item \strong{Group-level selection:} by determining whether \eqn{||\eta_{v}||_{2}=0}, we can know if the \eqn{v}th genetic variant has any effect at all.
#' \item \strong{Individual-level selection:} investigate whether the \eqn{v}th genetic variant has main effect, G\eqn{\times}E interaction or both, by determining which components in \eqn{\eta_{v}} has non-zero values.
#' }
#' If {structure="group"}, only group-level selection will be conducted on \eqn{||\eta_{v}||_{2}}; if {structure="individual"}, only individual-level selection will be conducted on each \eqn{\eta_{vu}}, (\eqn{u=1,\ldots,q}).
#'
#' This function also provides choices for the framework that is used.  If {func="QIF"}, variable selection will be conducted within the quadratic inference functions framework; if {func="GEE"}, variable selection will be
#' conducted within the generalized estimating equation framework.
#'
#' There are three options for the choice of the working correlation.  If {corr="exchangeable"}, the exchangeable working correlation will be applied; if {corr="AR-1"}, the AR-1 working correlation will be adopted; if {corr="independence"},
#' the independence working correlation will be used.
#' Please check the references for more details.
#' @examples
#' data("dat")
#' ##load the clinical covariates, environment factors, genetic factors and response from the
#' ##"dat" file
#' clin=dat$clin
#' if(is.null(clin)){t=0} else{t=dim(clin)[2]}
#' e=dat$e
#' u=dim(e)[2]
#' g=dat$g
#' y=dat$y
#' ##initial coefficient
#' beta0=dat$coef
#' ##true nonzero coefficients
#' index=dat$index
#' beta = springer(clin=clin, e, g, y,beta0,func="GEE",corr="independence",structure="bilevel",
#' lam1=dat$lam1, lam2=dat$lam2,maxits=30,tol=0.01)
#' ##only focus on the genetic main effects and gene-environment interactions
#' beta[1:(1+t+u)]=0
#' ##effects that have nonzero coefficients
#' pos = which(beta != 0)
#' ##true positive and false positive
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)

#'
#' @export

springer <- function(clin=NULL,e,g,y,beta0,func,corr,structure,lam1,lam2,maxits=30,tol=0.001){

  q=dim(e)[2]
  n=dim(y)[1]
  k=dim(y)[2]
  p1=dim(g)[2]

  e1=cbind(rep(1,dim(e)[1]),e)

  for (i in 1:p1) {
    e=cbind(e,g[,i]*e1)
  }

  if(is.null(clin)){
    t=0
    x=scale(e)}
  else{
    t=dim(clin)[2]
    x=scale(cbind(clin,e))
    }

  #==========================================reformat the data===============================#
  data=reformat(k,y,x)
  y=data$y
  x=data$x

  p=dim(x)[2]
  k=rep(k,n)

      converge=F
      iter=0
      beta.new=beta0

      if(func=="QIF"){func="Q"} else{func="G"}
      if(corr=="exchangeable"){corr="e"} else if(corr=="AR-1"){corr="a"} else{corr="i"}
      #=========================================PQIF======================================#
      while ((!converge) & (iter < maxits)) {
        beta = beta.new
        Score=ScoreU(n,k,y,x,p,beta,func,corr)
        U=Score$U
        dU=Score$dU

        E=penalty(x,n,t,p,q,beta,lam1,structure,p1,lam2)
        H=n*E
        matr1=dU + H
        matr2=U - H%*%beta

        beta.new=beta+NR(matr1,matr2)


        #mat[abs(mat)<0.000001]=0
        # save(mat,file = "mat.rda")
        # save.image(file = "example.RData")
        #beta.new = beta + ginv(matr1)%*%(U - H%*%beta)

        diff=mean(abs(beta.new-beta))
        converge = (diff < 1e-3)
        iter = iter+1
        #cat("iter",iter,"diff",diff,"\n")
      }
  coef=beta.new
  coef[abs(coef)<tol]=0
  return(coef)
}
