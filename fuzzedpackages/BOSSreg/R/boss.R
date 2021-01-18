#' Best Orthogonalized Subset Selection (BOSS).
#'
#'\itemize{
#'  \item Compute the solution path of BOSS and forward stepwise selection (FS).
#'  \item Compute various information criteria based on a heuristic degrees of freedom (hdf)
#'   that can serve as the selection rule to choose the subset given by BOSS.
#'   Only work when n>p.
#'}
#' @param x A matrix of predictors, with \code{nrow(x)=length(y)=n} observations and
#'   \code{ncol(x)=p} predictors. Intercept shall NOT be included.
#' @param y A vector of response variable, with \code{length(y)=n}.
#' @param intercept Logical, whether to include an intercept term. Default is TRUE.
#' @param hdf.ic.boss Logical, whether to calculate the heuristic degrees of freedom (hdf)
#'   and information criteria (IC) for BOSS. IC includes AIC, BIC, AICc, BICc, GCV,
#'   Cp. Note that if n<=p, \code{hdf.ic.boss=FALSE} no matter what. Default is TRUE.
#' @param mu True mean vector, used in the calculation of hdf. Default is NULL, and is estimated via
#'   least-squares (LS) regression of y upon x.
#' @param sigma True standard deviation of the error, used in the calculation of hdf. Default is NULL,
#'   and is estimated via least-squares (LS) regression of y upon x.
#' @param ... Extra parameters to allow flexibility. Currently none allows or requires, just for
#'   the convinience of call from other parent functions like cv.boss.
#'
#' @return
#' \itemize{
#'   \item beta_fs: A matrix of regression coefficients for all the subsets given by FS,
#'   from a null model until stop, with \code{nrow=p} and \code{ncol=min(n,p)+1}, where \code{min(n,p)} is
#'   the maximum number of steps performed.
#'   \item beta_boss: A matrix of regression coefficients for all the subsets given by
#'   BOSS, with \code{nrow=p} and \code{ncol=min(n,p)+1}. Note that unlike beta_fs and due to the nature of BOSS,
#'   the number of non-zero components in columns of beta_boss may not be unique, i.e.
#'   there maybe multiple columns corresponding to the same size of subset.
#'   \item steps_x: A vector of numbers representing which predictor joins at each step,
#'   with \code{length(steps)=min(n,p)}. The ordering is determined by the partial correlation between a predictor \eqn{x_j}
#'   and the response \code{y}.
#'   \item steps_q: A vector of numbers representing which predictor joins at each step in the orthogonal basis,
#'   with \code{length(steps)=min(n,p)}. BOSS takes the ordered predictors (ordering given in \code{steps_x}) and performs best
#'   subset regression upon their orthogonal basis, which is essentially ordering the orthogonalized predictors by their
#'   marginal correlations with the response \code{y}. For example, \code{steps_q=c(2,1)} indicates that the orthogonal basis of
#'   \code{x_2} joins first.
#'   \item hdf_boss: A vector of heuristic degrees of freedom (hdf) for BOSS, with
#'   \code{length(hdf_boss)=p+1}. Note that \code{hdf_boss=NULL} if n<=p or \code{hdf.ic.boss=FALSE}.
#'   \item IC_boss: A list of information criteria (IC) for BOSS, where each element
#'   in the list is a vector representing values of a given IC for each candidate subset
#'   of BOSS (or each column in beta_boss). The output IC includes AIC, BIC, AICc, BICc,
#'   GCV and Mallows' Cp. Note that each IC is calculated by plugging in hdf_boss.
#'
#' }
#'
#' @details This function computes the full solution path given by BOSS and FS on a given
#'   dataset (x,y) with n observations and p predictors. In the case where n>p, it also calculates
#'   the heuristic degrees of freedom for BOSS, and various information criteria, which can further
#'   be used to select the subset from the candidates. Please refer to the Vignette
#'   for implementation details and Tian et al. (2019) for methodology details (links are given below).
#'
#' @author Sen Tian
#' @references
#' \itemize{
#'   \item Tian, S., Hurvich, C. and Simonoff, J. (2019), On the Use of Information Criteria
#'   for Subset Selection in Least Squares Regression. https://arxiv.org/abs/1911.10191
#'   \item BOSSreg Vignette https://github.com/sentian/BOSSreg/blob/master/r-package/vignettes/BOSSreg.pdf
#' }
#' @seealso \code{predict} and  \code{coef} methods for "boss" object, and the \code{cv.boss} function
#' @example R/example/eg.boss.R
#' @useDynLib BOSSreg
#' @importFrom Rcpp sourceCpp
#' @export
boss <- function(x, y, intercept=TRUE, hdf.ic.boss=TRUE, mu=NULL, sigma=NULL, ...){
  n = dim(x)[1]
  p = dim(x)[2]

  maxstep = min(n, p)
  varnames = colnames(x)
  # standardize x (mean 0 and norm 1) and y (mean 0)
  std_result = std(x, y, intercept)
  x = std_result$x_std
  y = std_result$y_std
  mean_x = std_result$mean_x
  mean_y = std_result$mean_y
  sd_demanedx = std_result$sd_demeanedx

  guideQR_result = guideQR(x, y, maxstep)
  Q = guideQR_result$Q
  R = guideQR_result$R
  steps_x = as.numeric(guideQR_result$steps)

  # coefficients
  z = t(Q) %*% y

  # fs
  beta_q = matrix(rep(z, maxstep), nrow=maxstep, byrow=F)
  beta_q = beta_q * upper.tri(beta_q, diag=TRUE)
  beta_q = cbind(0, beta_q)
  if(n<p){
    steps_expand = c(steps_x, setdiff(seq(1, p),steps_x))
    beta_fs = rbind(backsolve(R, beta_q)[order(steps_x),], matrix(0, nrow=p-n, ncol=maxstep+1))[order(steps_expand), ]
  }else{
    beta_fs = backsolve(R, beta_q)[order(steps_x), ]
  }
  # scale back
  beta_fs = diag(1/sd_demanedx) %*% beta_fs
  #beta_fs = cbind(0,beta_fs)
  if(intercept){
    beta_fs = rbind((mean_y - mean_x %*% beta_fs), beta_fs)
  }
  beta_fs = Matrix::Matrix(beta_fs, sparse=TRUE)

  beta_boss = hdf_result = IC_result = NULL


  # boss
  order_q = order(-z^2)
  steps_q = steps_x[order_q]
  row_i = rep(order_q, times=seq(p,1))
  col_j = unlist(lapply(2:(p+1), function(xx){seq(xx,p+1)}))
  beta_q = Matrix::sparseMatrix(row_i, col_j, x=z[row_i], dims=c(p, p+1))

  # old version
  # beta_q = matrix(0, nrow=p, ncol=maxstep+1)
  # order_q = order(-z^2)
  # steps_q = steps_x[order_q]
  # for(j in 1:maxstep){
  #   beta_q[order_q[1:j], (j+1)] = z[order_q[1:j]]
  # }

  # project back and change the order
  if(n<p){
    beta_boss = rbind(backsolve(R, beta_q), matrix(0, nrow=p-n, ncol=maxstep+1))[order(steps_expand), ]
  }else{
    beta_boss = backsolve(R, beta_q)[order(steps_x), ]
  }
  # scale back
  beta_boss = diag(1/sd_demanedx) %*% beta_boss
  if(intercept){
    beta_boss = rbind((mean_y - mean_x %*% beta_boss), beta_boss)
  }
  # beta_boss = Matrix::Matrix(beta_boss, sparse=TRUE)

  # hdf and IC
  if(n<=p & hdf.ic.boss){
    warning('hdf not available when n<=p')
  }else if(hdf.ic.boss){
    hdf_result = calc.hdf(Q, y, sigma, mu)
    if(intercept){
      hdf_result$hdf = hdf_result$hdf + 1
    }
    if(is.null(sigma)){
      IC_result = calc.ic.all(beta_q, Q, y, hdf_result$hdf, hdf_result$sigma)
    }else{
      IC_result = calc.ic.all(beta_q, Q, y, hdf_result$hdf, sigma)
    }
  }

  # take care the variable names
  if(is.null(varnames)){
    if(intercept){
      rownames(beta_fs) = rownames(beta_boss) = c('intercept', paste('X',seq(1,p),sep=''))
    }else{
      rownames(beta_fs) = rownames(beta_boss) = paste('X',seq(1,p),sep='')
    }
    names(steps_x) = paste('X',steps_x,sep='')
    names(steps_q) = paste('X',steps_q,sep='')
  }else{
    if(intercept){
      rownames(beta_fs) = rownames(beta_boss) = c('intercept', varnames)
    }else{
      rownames(beta_fs) = rownames(beta_boss) = varnames
    }
    names(steps_x) = varnames[steps_x]
    names(steps_q) = varnames[steps_q]
  }
  # output
  out = list(beta_fs=beta_fs,
             beta_boss=beta_boss,
             steps_x=steps_x,
             steps_q=steps_q,
             hdf_boss=hdf_result$hdf,
             IC_boss=IC_result,
             call=list(intercept=intercept))
  class(out) = 'boss'
  invisible(out)
}





#' Select coefficient vector(s) for BOSS.
#'
#' This function returns the optimal coefficient vector of BOSS selected by AICc
#' (by default) or other types of information criterion.
#'
#' @param object The boss object, returned from calling the \code{boss} function.
#' @param ic Which information criterion is used to select the optimal coefficient vector for BOSS.
#' The default is AICc-hdf.
#' @param select.boss The index (or indicies) of columns in the coefficient matrix for which
#' one wants to select. By default (NULL) it's selected by the information criterion specified in
#' 'ic'.
#' @param ... Extra arguments (unused for now)
#'
#' @return The chosen coefficient vector(s) for BOSS.
#'
#' @details If \code{select.boss} is specified, the function returns
#' corresponding column(s) in the coefficient matrix.
#'
#' If \code{select.boss} is unspecified, the function returns the optimal coefficient
#' vector selected by AICc-hdf (other choice of IC can be specified in the argument \code{ic}).
#' The only exception is when n>=p, where hdf is not well defined, and the entire coefficient matrix
#' is returned.
#'
#' @example R/example/eg.boss.R
#' @importFrom stats coef
#' @export
coef.boss <- function(object, ic=c('aicc','bicc','aic','bic','gcv','cp'), select.boss=NULL, ...){
  # # for fs, return the full coef matrix if not specified the columns
  # if(is.null(select.fs)){
  #   select.fs = 1:ncol(object$beta_fs)
  # }else if(select.fs == 0){
  #   select.fs = 1:ncol(object$beta_fs)
  # }
  # select.fs[select.fs > ncol(object$beta_fs)] = ncol(object$beta_fs)
  # beta_fs_opt = object$beta_fs[, select.fs, drop=FALSE]

  # for boss, the default is to return coef selected by AICc
  if(is.null(select.boss)){
    if(is.null(object$IC_boss)){
      # if we are in the case where n>=p
      if(dim(object$beta_boss)[1]+1 >= dim(object$beta_boss)[2]){
        warning('hdf does not work when n>=p, full coef matrix returned')
      }else{
        # this is where hdf.ic.boss is flagged FALSE
        warning("rerun boss with argument 'hdf.ic.boss=TRUE', and call coef.boss again")
      }
      select.boss = 1:ncol(object$beta_boss)
    }else{
      ic = match.arg(ic)
      select.boss = which.min(object$IC_boss[[ic]])
    }
  }else if(select.boss == 0){
    select.boss = 1:ncol(object$select.boss)
  }
  select.boss[select.boss > ncol(object$beta_boss)] = ncol(object$beta_boss)
  beta_boss_opt = object$beta_boss[, select.boss, drop=FALSE]

  # return(list(fs=beta_fs_opt, boss=beta_boss_opt))
  return(beta_boss_opt)
}

#' Prediction given new data entries.
#'
#' This function returns the prediction(s) given new observation(s), for BOSS,
#' where the optimal coefficient vector is chosen via certain selection rule.
#'
#' @param object The boss object, returned from calling 'boss' function.
#' @param newx A new data entry or several entries. It can be a vector, or a matrix with
#' \code{nrow(newx)} being the number of new entries and \code{ncol(newx)=p} being the
#' number of predictors. The function takes care of the intercept, NO need to add \code{1}
#' to \code{newx}.
#' @param ... Extra arguments to be plugged into \code{coef}, such as \code{select.boss},
#' see the description of \code{coef.boss} for more details.
#'
#' @return The prediction(s) for BOSS.
#'
#' @details The function basically calculates \eqn{x * coef}, where \code{coef}
#' is a coefficient vector chosen by a selection rule. See more details about the default
#' and available choices of the selection rule in the description of \code{coef.boss}.
#'
#' @example R/example/eg.boss.R
#' @importFrom stats predict
#' @export
predict.boss <- function(object, newx, ...){
  # coefficients
  # coef_result = coef(object, ...)
  # beta_fs_opt = coef_result$fs
  # beta_boss_opt = coef_result$boss
  beta_boss_opt = coef(object, ...)

  # make newx a matrix
  # if newx is an array or a column vector, make it a row vector
  if(is.null(dim(newx))){
    newx = matrix(newx, nrow=1)
  }else if(dim(newx)[2] == 1){
    newx = t(newx)
  }
  # if intercept, add 1 to newx
  if(object$call$intercept){
    newx = cbind(rep(1,nrow(newx)), newx)
  }

  # check the dimension
  # if(ncol(newx) != nrow(beta_fs_opt)){
  #   stop('Mismatch dimension of newx and coef for FS. Note do NOT add 1 to newx when intercept=TRUE')
  # }else{
  #   mu_fs_opt = newx %*% beta_fs_opt
  # }
  if(is.null(beta_boss_opt)){
    mu_boss_opt = NULL
  }else{
    if(ncol(newx) != nrow(beta_boss_opt)){
      stop('Mismatch dimension of newx and coef for BOSS. Note do NOT add 1 to newx when intercept=TRUE')
    }else{
      mu_boss_opt = newx %*% beta_boss_opt
    }
  }
  # return(list(fs=mu_fs_opt, boss=mu_boss_opt))
  return(mu_boss_opt)
}








