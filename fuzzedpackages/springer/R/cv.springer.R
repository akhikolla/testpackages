#' @useDynLib springer, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' k-folds cross-validation for springer
#'
#' This function conducts k-fold cross-validation for springer and returns the optimal values of the tuning parameters.
#' @param clin a matrix of clinical covariates. The default value is NULL. Whether to include the clinical covariates is decided by user.
#' @param e a matrix of environment factors.
#' @param g a matrix of genetic factors.
#' @param y the longitudinal response.
#' @param beta0 the initial value for the coefficient vector.
#' @param lambda1 a user-supplied sequence of \eqn{\lambda_{1}} values, which serves as a tuning parameter for the individual-level penalty.
#' @param lambda2 a user-supplied sequence of \eqn{\lambda_{2}} values, which serves as a tuning parameter for the group-level penalty.
#' @param nfolds the number of folds for cross-validation.
#' @param func the framework to obtain the score equation.  Two choices are available: "GEE" and "QIF".
#' @param corr the working correlation structure adopted in the estimation algorithm. The springer provides three choices for the
#' working correlation structure: exchangeable, AR-1,and independence.
#' @param structure Three choices are available for structured variable selection. "bilevel" for sparse-group selection on both group-level and individual-level. "group" for selection on group-level only. "individual" for selection on individual-level only.
#' @param maxits the maximum number of iterations that is used in the estimation algorithm. The default value is 30.
#' @param tol The tolerance level. Coefficients with absolute values that are smaller than the tolerance level will be set to zero. The adhoc value can be chosen as 0.001.
#' @details
#' For bi-level sparse group selection, cv.springer returns two optimal tuning parameters,
#' \eqn{\lambda_{1}} and \eqn{\lambda_{2}}; for group-level selection, this function returns the optimal \eqn{\lambda_{2}} with \eqn{\lambda_{1}}=0;
#' for individual-level selection, this function returns the optimal \eqn{\lambda_{1}} with \eqn{\lambda_{2}}=0.
#' @return an object of class "cv.springer" is returned, with is a list with components below:
#' \item{lam1}{the optimal \eqn{\lambda_{1}}.}
#' \item{lam2}{the optimal \eqn{\lambda_{2}}.}
#'
#' @export

cv.springer <- function(clin=NULL, e, g, y, beta0, lambda1, lambda2, nfolds, func, corr, structure, maxits=30, tol=0.001){
  n=dim(y)[1]
  q=dim(e)[2]
  p1=dim(g)[2]
  len1=length(lambda1)
  len2=length(lambda2)

  foldsize=round(n/nfolds,0)
  folds=rep(1,n-foldsize*(n-1))
  for (i in 2:nfolds) {
    folds=c(folds,rep(i,n/nfolds))
  }
  folds=sample(folds)
  pred=matrix(0,len1,len2)
  for (i in 1:len1) {
    lam1=lambda1[i]
    for (j in 1:len2){
      lam2=lambda2[j]
      mse=0
      for (cv in 1:nfolds) {
        #Segement your data by fold using the which() function
        testIndexes <- which(folds==cv,arr.ind=TRUE)

        e.test=e[testIndexes,]
        g.test=g[testIndexes,]
        y.test=y[testIndexes,]

        e.train=e[-testIndexes,]
        g.train=g[-testIndexes,]
        y.train=y[-testIndexes,]
        e.train=as.matrix(e.train)
        g.train=as.matrix(g.train)
        y.train=as.matrix(y.train)

        if(is.null(clin)){
          clin.train=clin.test=NULL
        }
        else{
          clin.test=clin[testIndexes,]
          clin.train=clin[-testIndexes,]
        }
        beta=springer(clin=clin.train,e.train, g.train, y.train, beta0,func,corr,structure,lam1,lam2,maxits,tol)

        e1=cbind(rep(1,dim(e.test)[1]),e.test)

        for (i in 1:p1) {
          e.test=cbind(e.test,g.test[,i]*e1)
        }

        if(is.null(clin.test)){
          x.test=scale(e.test)
        }
        else{
          x.test=scale(cbind(clin.test,e.test))
        }

        data.test=reformat(y.test, x.test)
        x.test=data.test$x
        y.test=data.test$y
        mu=x.test%*%beta
        mse=mse+mean((y.test-mu)^2)
      }
      pred[i,j]=mse/nfolds
    }
  }
  lamb1=lambda1[which(pred==min(pred),arr.ind = TRUE)[1]]
  lamb2=lambda2[which(pred==min(pred),arr.ind = TRUE)[2]]
  return(list("lam1"=lamb1,"lam2"=lamb2))
}
