#' Cross-validation for Best Orthogonalized Subset Selection (BOSS) and Forward Stepwise Selection (FS).
#'
#' @param x A matrix of predictors, see \code{boss}.
#' @param y A vector of response variable, see \code{boss}.
#' @param n.folds The number of cross validation folds. Default is 10.
#' @param n.rep The number of replications of cross validation. Default is 1.
#' @param intercept Logical, whether to fit an intercept term. Default is TRUE.
#' @param ... Arguments to \code{boss}, such as \code{hdf.ic.boss}.
#'
#' @return
#' \itemize{
#'   \item boss: An object \code{boss} that fits on the full dataset.
#'   \item n.folds: The number of cross validation folds.
#'   \item cvm.fs: Mean OOS deviance for each candidate given by FS.
#'   \item cvm.boss: Mean OSS deviance for each candidate given by BOSS.
#'   \item i.min.fs: The index of minimum cvm.fs.
#'   \item i.min.boss: The index of minimum cvm.boss.
#' }
#'
#' @details This function fits BOSS and FS (\code{boss}) on the full dataset, and performs \code{n.folds}
#'   cross-validation. The cross-validation process can be repeated \code{n.rep} times to evaluate the
#'   out-of-sample (OOS) performance for the candidate subsets given by both methods.
#'
#' @author Sen Tian
#' @references
#' \itemize{
#'   \item Tian, S., Hurvich, C. and Simonoff, J. (2019), On the Use of Information Criteria
#'   for Subset Selection in Least Squares Regression. https://arxiv.org/abs/1911.10191
#'   \item BOSSreg Vignette https://github.com/sentian/BOSSreg/blob/master/r-package/vignettes/BOSSreg.pdf
#' }
#' @seealso \code{predict} and \code{coef} methods for \code{cv.boss} object, and the \code{boss} function
#' @example R/example/eg.cv.boss.R
#' @export
cv.boss <- function(x, y, n.folds=10, n.rep=1, intercept=TRUE, ...){
  # # arguments
  argu = list(...)
  # argu_boss = c('intercept', 'hdf.ic.boss') # arguments that boss accepts
  # # arguments that user specify but unused
  # argu_unused = setdiff(names(argu), argu_boss)
  # if(length(argu_unused) > 0){
  #   warning(paste(argu_unused, ' are not valid arguments for boss, check spelling maybe?', sep=''))
  # }

  # overide hdf.ic.boss option in '...', to be used in CV
  # boss.nohdf <- function(x, y, intercept, hdf.ic.boss) boss(x, y, intercept, hdf.ic.boss=FALSE)

  # start the CV process
  n = dim(x)[1]
  p = dim(x)[2]
  maxstep = trunc(min(n - n/n.folds, p))
  if(maxstep < p){
    warning('the number of observations in each fold, does not allow evaluating the full path, some large sets of variables are ignored')
  }

  # matrix to store the CV error
  cv_rep_boss = cv_rep_fs = matrix(NA, nrow=n.rep, ncol=maxstep+1)

  for(replication in 1:n.rep){
    fold.index = sample(rep(1:n.folds, length.out=n)) # randomly assign a fold to each observation
    cv_tmp_boss = cv_tmp_fs = matrix(NA, nrow=n.folds, ncol=maxstep+1)

    for(fold in 1:n.folds){
      # split the training and testing sets
      test.index = which(fold.index==fold)
      x.test = x[test.index, , drop=FALSE]
      y.test = y[test.index]
      x.train = x[-test.index, , drop=FALSE]
      y.train = y[-test.index]
      boss_result = boss(x.train, y.train, intercept, hdf.ic.boss=FALSE)
      beta_fs = boss_result$beta_fs
      beta_boss = boss_result$beta_boss
      # if intercept
      if(intercept){
        x.test = cbind(rep(1,nrow(x.test)), x.test)
      }
      cv_tmp_fs[fold, ] = Matrix::colMeans(sweep(x.test%*%beta_fs, 1, y.test, '-')^2)
      cv_tmp_boss[fold, ] = Matrix::colMeans(sweep(x.test%*%beta_boss, 1, y.test, '-')^2)

    }
    cv_rep_fs[replication, ] = Matrix::colMeans(cv_tmp_fs)
    cv_rep_boss[replication, ] = Matrix::colMeans(cv_tmp_boss)
  }

  cv_fs = Matrix::colMeans(cv_rep_fs)
  cv_boss = Matrix::colMeans(cv_rep_boss)


  # fit on the full sample
  boss_result = boss(x, y, intercept, ...)

  # output
  out = list(boss=boss_result,
             n.folds=n.folds,
             cvm.fs=cv_fs,
             cvm.boss=cv_boss,
             i.min.fs=which.min(cv_fs),
             i.min.boss=which.min(cv_boss),
             call=list(intercept=intercept))
  class(out) = 'cv.boss'
  invisible(out)
}


#' Select coefficient vector based on cross-validation for BOSS or FS.
#'
#' This function returns coefficient vector that minimizes out-of-sample (OOS) cross
#' validation score.
#'
#' @param object The cv.boss object, returned from calling \code{cv.boss} function.
#' @param method It can either be 'fs' or 'boss'. The default is 'boss'.
#' @param ... Extra arguments (unused for now).
#'
#' @return The chosen coefficient vector for BOSS or FS.
#'
#' @example R/example/eg.cv.boss.R
#' @importFrom stats coef
#' @export
coef.cv.boss <- function(object, method=c('boss', 'fs'), ...){
  # coef_result = coef(object$boss, select.fs=object$i.min.fs, select.boss=object$i.min.boss)
  # beta_fs_opt = coef_result$fs
  # beta_boss_opt = coef_result$boss
  # return(list(fs=beta_fs_opt, boss=beta_boss_opt))

  if(match.arg(method) == 'fs'){
    beta_fs_opt = object$boss$beta_fs[, object$i.min.fs, drop=FALSE]
    return(beta_fs_opt)
  }else{
    beta_boss_opt = coef(object$boss, select.boss=object$i.min.boss)
    return(beta_boss_opt)
  }
}

#' Prediction given new data entries.
#'
#' This function returns the prediction(s) given new observation(s) for BOSS or FS,
#' where the optimal coefficient vector is chosen via cross-validation.
#'
#' @param object The cv.boss object, returned from calling \code{cv.boss} function.
#' @param newx A new data entry or several entries. It can be a vector, or a matrix with
#' \code{nrow(newx)} being the number of new entries and \code{ncol(newx)=p} being the
#' number of predictors. The function takes care of the intercept, NO need to add \code{1}
#' to \code{newx}.
#' @param ... Extra arguments to be plugged into \code{coef}, such as \code{method},
#' see the description of \code{coef.cv.boss} for more details.
#'
#' @return The prediction for BOSS or FS.
#'
#' @example R/example/eg.cv.boss.R
#' @importFrom stats predict
#' @export
predict.cv.boss <- function(object, newx, ...){
  # predict_result = predict(object$boss, newx, select.fs=object$i.min.fs, select.boss=object$i.min.boss)
  # mu_fs_opt = predict_result$fs
  # mu_boss_opt = predict_result$boss
  # return(list(fs=mu_fs_opt, boss=mu_boss_opt))

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

  beta_opt = coef(object, ...)

  # check the dimension
  if(ncol(newx) != nrow(beta_opt)){
    stop('Mismatch dimension of newx and coef. Note do NOT add 1 to newx when intercept=TRUE')
  }else{
    mu_opt = newx %*% beta_opt
  }
  return(mu_opt)

}
