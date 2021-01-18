#' Jittering with Over/Under Sampling
#'
#' Perform probability estimation using jittering with over or undersampling.
#'
#' @param X A matrix of continuous predictors.
#' @param y A vector of responses with entries in \code{c(-1, 1)}.
#' @param delta An integer (greater than 3) to control the number of quantiles to
#'        estimate:
#' @param class_func Function to perform classification.  This function definition must be
#'        exactly of the form \code{class_func(X, y)} where X is a matrix and y is a
#'        vector with entries in \code{c(-1, 1)}, and it must return an object on which
#'        \code{pred_func} can create predictions.  See examples.
#' @param pred_func Function to create predictions.  This function definition must be
#'        exactly of the form \code{pred_func(fit_obj, X)} where \code{fit_obj}
#'        is an object returned by class_func and X is a matrix of new data
#'        values, and it must return a vector with entries in \code{c(-1, 1)}.
#'        See examples.
#' @param type Type of sampling: "over" for oversampling,  or "under" for
#'        undersampling.
#' @param nu The amount of noise to apply to predictors when oversampling data.
#'        The noise level is controlled by \code{nu * sd(X[,j])} for each
#'        predictor - the default of \code{nu = 1} works well.  Such "jittering"
#'        of the predictors is essential when applying \code{jous} to boosting
#'        type methods.
#' @param X_pred A matrix of predictors for which to form probability estimates.
#' @param keep_models Whether to store all of the models used to create
#'        the probability estimates.  If \code{type=FALSE}, the user will need
#'        to re-run \code{jous} when creating probability estimates for test data.
#' @param verbose If \code{TRUE}, print the function's progress to the terminal.
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit models.
#'        Must register parallel before hand, such as \code{doParallel}.  See
#'        examples below.
#' @param packages If \code{parallel = TRUE}, a vector of strings containing the
#'        names of any packages used in \code{class_func} or \code{pred_func}.
#'        See examples below.
#' @return Returns a list containing information about the
#' parameters used in the \code{jous} function call, as well as the following
#' additional components:
#' \item{q}{The vector of target quantiles estimated by \code{jous}.  Note that
#' the estimated probabilities will be located at the midpoints of the values in
#' \code{q}.}
#' \item{phat_train}{The in-sample probability estimates \eqn{p(y=1|x)}.}
#' \item{phat_test}{Probability estimates for the optional test data in \code{X_test}}
#' \item{models}{If \code{keep_models=TRUE}, a list of models fitted to
#' the resampled data sets.}
#' \item{confusion_matrix}{A confusion matrix for the in-sample fits.}
#'
#' @note The \code{jous} function runs the classifier \code{class_func} a total
#' of \code{delta} times on the data, which can be computationally expensive.
#' Also,\code{jous} cannot yet be applied to categorical predictors - in the
#' oversampling case, it is not clear how to "jitter" a categorical variable.
#'
#' @references Mease, D., Wyner, A. and Buja, A. (2007). Costweighted
#' boosting with jittering and over/under-sampling:
#' JOUS-boost. J. Machine Learning Research 8 409-439.
#'
#' @examples
#' \dontrun{
#' # Generate data from Friedman model #
#' set.seed(111)
#' dat = friedman_data(n = 500, gamma = 0.5)
#' train_index = sample(1:500, 400)
#'
#' # Apply jous to adaboost classifier
#' class_func = function(X, y) adaboost(X, y, tree_depth = 2, n_rounds = 200)
#' pred_func = function(fit_obj, X_test) predict(fit_obj, X_test)
#'
#' jous_fit = jous(dat$X[train_index,], dat$y[train_index], class_func,
#'                 pred_func, keep_models = TRUE)
#' # get probability
#' phat_jous = predict(jous_fit, dat$X[-train_index, ], type = "prob")
#'
#' # compare with probability from AdaBoost
#' ada = adaboost(dat$X[train_index,], dat$y[train_index], tree_depth = 2,
#'                n_rounds = 200)
#' phat_ada = predict(ada, dat$X[train_index,], type = "prob")
#'
#' mean((phat_jous - dat$p[-train_index])^2)
#' mean((phat_ada - dat$p[-train_index])^2)
#'
#' ## Example using parallel option
#'
#' library(doParallel)
#' cl <- makeCluster(4)
#' registerDoParallel(cl)
#'
#' # n.b. the packages='rpart' is not really needed here since it gets
#' # exported automatically by JOUSBoost, but for illustration
#' jous_fit = jous(dat$X[train_index,], dat$y[train_index], class_func,
#'                 pred_func, keep_models = TRUE, parallel = TRUE,
#'                 packages = 'rpart')
#' phat = predict(jous_fit, dat$X[-train_index,], type = 'prob')
#' stopCluster(cl)
#'
#' ## Example using SVM
#'
#' library(kernlab)
#' class_func = function(X, y) ksvm(X, as.factor(y), kernel = 'rbfdot')
#' pred_func = function(obj, X) as.numeric(as.character(predict(obj, X)))
#' jous_obj = jous(dat$X[train_index,], dat$y[train_index], class_func = class_func,
#'            pred_func = pred_func, keep_models = TRUE)
#' jous_pred = predict(jous_obj, dat$X[-train_index,], type = 'prob')
#' }
#' @export
jous = function(X, y,
                class_func,
                pred_func,
                type=c("under", "over"),
                delta=10,
                nu=1,
                X_pred=NULL,
                keep_models=FALSE,
                verbose=FALSE,
                parallel=FALSE,
                packages=NULL){

  # check data types
  if(!all(y %in% c(-1,1)))
    stop("y must take values in -1, 1")

  if(!is.matrix(X))
    stop("X must be a matrix")

  if(delta < 3)
    stop("delta must be an integer greater than 2")

  # extract tpye of sampling
  type = match.arg(type)

  # warn if no packages specified and parallel = T
  if(parallel & is.null(packages)){
    message(paste0('"parallel" = TRUE, but no packages specified for export.',
                   '  If function fails, specify a non-null value for ',
                   '"packages."'))
    packages=''
  }

  ix_pos = which(y == 1)
  ix_neg = which(y == -1)
  ncuts = delta - 1

  q = (1:ncuts)/delta
  # if the 0.5 isn't in q, insert it
  if(!any(q == 0.5)){
    q <- c(q[q < 0.5], 0.5, q[q > 0.5])
    ncuts = ncuts + 1
  }

  ## Fit models over tilted data ##
  models = list()
  if(type == 'over'){
    ix = index_over(ix_pos, ix_neg, q)
    col_stds = apply(X, 2, stats::sd, na.rm=TRUE)
    X_ = sapply(1:ncol(X), function(i) X[,i] +
                  stats::runif(n=nrow(X), -col_stds[i]*nu, col_stds[i]*nu))
  } else{
    ix = index_under(ix_pos, ix_neg, q, delta)
    X_ = X # redundant copy, but will simplify code a bit
  }

  # loop over data-sets
  if(parallel){
    i = NULL
    models = foreach::`%dopar%`(foreach::foreach(i = seq(ncuts), .inorder=T,
                              .packages=packages),
                              {
                                ix_temp = c(ix$ix_neg_cut[[i]],
                                            ix$ix_pos_cut[[i]])
                                if(type == 'over'){
                                  class_func(rbind(X, X_[ix_temp,]),
                                             c(y, y[ix_temp]))
                                  } else{
                                    class_func(X[ix_temp,], y[ix_temp])
                                  }
                                })
  }else{
    for(i in seq(ncuts)){
      ix_temp = c(ix$ix_neg_cut[[i]], ix$ix_pos_cut[[i]])
      if(type == 'over'){
        models[[i]] = class_func(rbind(X, X_[ix_temp,]), c(y, y[ix_temp]))
      } else{
        models[[i]] = class_func(X[ix_temp,], y[ix_temp])
      }
      if(verbose) cat('Done with iteration ', i, ' of ', ncuts, '\n')
    }
  }

  # create jous object
  jous_obj = list(delta=delta, nu=nu, q=q, type=type, models=models,
                  pred_func = pred_func, packages=packages, parallel=parallel)
  class(jous_obj) = "jous"

  # in sample
  jous_obj$phat_train = stats::predict(jous_obj, X, type="prob")

  # create confusion matrix for in-sample fits
  yhat = ifelse(jous_obj$phat_train < 0.5, -1, 1)
  jous_obj$confusion_matrix = table(y, yhat)

  # out of sample
  if(!is.null(X_pred))
    jous_obj$phat_test = stats::predict(jous_obj, X_pred, type="prob")

  if(!keep_models)
    jous_obj$models = NULL

  jous_obj

}

#' Create predictions
#'
#' Makes a prediction on new data for a given fitted \code{jous} model.
#' @param object An object of class \code{jous} returned by the \code{jous} function.
#' @param X A design matrix of predictors.
#' @param type The type of prediction to return.  If \code{type="response"}, a
#'        class label of -1 or 1 is returned.  If \code{type="prob"}, the
#'        probability \eqn{p(y=1|x)} is returned.
#' @param ... \dots
#'
#' @return Returns a vector of class predictions if \code{type="response"}, or a
#'          vector of class probabilities \eqn{p(y=1|x)} if \code{type="prob"}.
#'
#' @examples
#' \dontrun{
#' # Generate data from Friedman model #
#' set.seed(111)
#' dat = friedman_data(n = 500, gamma = 0.5)
#' train_index = sample(1:500, 400)
#'
#' # Apply jous to adaboost classifier
#' class_func = function(X, y) adaboost(X, y, tree_depth = 2, n_rounds = 100)
#' pred_func = function(fit_obj, X_test) predict(fit_obj, X_test)
#'
#' jous_fit = jous(dat$X[train_index,], dat$y[train_index], class_func,
#'                 pred_func, keep_models=TRUE)
#' # get class prediction
#' yhat = predict(jous_fit, dat$X[-train_index, ])
#' # get probability estimate
#' phat = predict(jous_fit, dat$X[-train_index, ], type="prob")
#' }
#' @export
predict.jous = function(object, X, type=c("response", "prob"), ...){

  # handle args
  type = match.arg(type)
  if(is.null(object$models))
    stop("No saved models in your jous object.  Rerun with keep_models = TRUE")

  delta = object$delta
  q = object$q

  ## calculate predictions for each classifier
  if(object$parallel){
    i = NULL
    pred_mat = foreach::`%dopar%`(foreach::foreach(i = seq_along(object$models),
                                                   .inorder=T,
                                                   .packages = object$packages,
                                                   .combine=cbind),
                                  {
                                    object$pred_func(object$models[[i]], X)
                                  })
  } else{
    pred_mat = sapply(object$models, function(z) object$pred_func(z, X))
  }

  if(!all(pred_mat %in% c(-1,1)))
    stop("Your prediction function must return values only in -1, 1")

  median_loc = which(q == 0.5) - 1 # 0 based indexing here
  phat = grid_probs(pred_mat, q, delta, median_loc)

  # handle response type
  if(type == "response"){
    2*(phat > 0.5) - 1
  } else if(type =="prob"){
    phat
  }

}

#' Print a summary of \code{jous} fit.
#' @param x A \code{jous} object.
#' @param ... \dots
#' @return Printed summary of the fit
#' @export
print.jous = function(x, ...){
  cat('jous fit: \n')
  cat('    type: ', x$type, '\n')
  cat('    delta: ', x$delta, '\n')
  cat('\n In-sample confusion matrix:\n')
  print(x$confusion_matrix)
}

